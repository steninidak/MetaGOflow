#!/usr/bin/env python3
"""Create an RO-Crate 1.1 ZIP for a completed MetaGOflow Nextflow run."""

import argparse
import hashlib
import json
import mimetypes
from pathlib import Path
from datetime import datetime, timezone
from zipfile import ZIP_DEFLATED, ZIP_STORED, ZipFile


def checksum(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_entity(identifier, source):
    media_type = mimetypes.guess_type(source.name)[0] or "application/octet-stream"
    return {
        "@id": identifier,
        "@type": "File",
        "name": source.name,
        "encodingFormat": media_type,
        "contentSize": source.stat().st_size,
        "sha256": checksum(source),
    }


def iter_payload(root, excluded):
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.resolve() not in excluded:
            yield path


def add_file(archive, source, destination):
    # Already-compressed bioinformatics and archive formats should not be
    # recompressed; doing so wastes substantial CPU with no useful reduction.
    stored_suffixes = {".gz", ".zip", ".jar", ".sif", ".bam", ".cram"}
    compression = ZIP_STORED if source.suffix.lower() in stored_suffixes else ZIP_DEFLATED
    archive.write(source, destination, compress_type=compression)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results", required=True, type=Path)
    parser.add_argument("--workflow", required=True, type=Path)
    parser.add_argument("--params", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--run-name", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--nextflow-version", required=True)
    args = parser.parse_args()

    results = args.results.resolve()
    workflow = args.workflow.resolve()
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)

    graph = [
        {
            "@id": "ro-crate-metadata.json",
            "@type": "CreativeWork",
            "about": {"@id": "./"},
            "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"},
        },
        {
            "@id": "./",
            "@type": "Dataset",
            "name": f"MetaGOflow Nextflow results: {args.run_name}",
            "datePublished": datetime.now(timezone.utc).isoformat(),
            "license": {"@id": "https://creativecommons.org/licenses/by/4.0/"},
            "publisher": {"@id": "https://ror.org/0038zss60"},
            "mentions": [{"@id": "#workflow-run"}],
            "hasPart": [],
        },
        {
            "@id": "https://ror.org/0038zss60",
            "@type": "Organization",
            "name": "European Marine Biological Resource Centre",
            "url": "https://www.embrc.eu/",
        },
        {
            "@id": "https://creativecommons.org/licenses/by/4.0/",
            "@type": "CreativeWork",
            "name": "Creative Commons Attribution 4.0 International",
        },
        {
            "@id": "#workflow-run",
            "@type": "CreateAction",
            "name": f"MetaGOflow Nextflow execution {args.run_name}",
            "actionStatus": {"@id": "http://schema.org/CompletedActionStatus"},
            "instrument": {"@id": "workflow/main.nf"},
            "identifier": args.run_id,
            "result": {"@id": "results/"},
            "object": {"@id": "workflow/params.json"},
            "softwareRequirements": f"Nextflow {args.nextflow_version}",
        },
        {
            "@id": "results/",
            "@type": "Dataset",
            "name": "MetaGOflow analysis results",
            "hasPart": [],
        },
        {
            "@id": "workflow/",
            "@type": "Dataset",
            "name": "MetaGOflow Nextflow workflow sources",
            "hasPart": [],
        },
    ]
    root = graph[1]
    results_dataset = graph[5]
    workflow_dataset = graph[6]
    root["hasPart"].extend([{"@id": "results/"}, {"@id": "workflow/"}])

    metadata = {
        "@context": "https://w3id.org/ro/crate/1.1/context",
        "@graph": graph,
    }

    excluded = {output}
    with ZipFile(output, "w", allowZip64=True) as archive:
        for source in iter_payload(results, excluded):
            relative = source.relative_to(results)
            identifier = f"results/{relative.as_posix()}"
            add_file(archive, source, identifier)
            graph.append(file_entity(identifier, source))
            results_dataset["hasPart"].append({"@id": identifier})

        workflow_files = [
            workflow / "main.nf",
            workflow / "nextflow.config",
            workflow / "nextflow_schema.json",
        ]
        workflow_files.extend(sorted((workflow / "nextflow").rglob("*.nf")))
        for source in workflow_files:
            if not source.is_file():
                continue
            relative = source.relative_to(workflow)
            identifier = f"workflow/{relative.as_posix()}"
            add_file(archive, source, identifier)
            entity = file_entity(identifier, source)
            if relative.as_posix() == "main.nf":
                entity["@type"] = ["File", "SoftwareSourceCode", "ComputationalWorkflow"]
                entity["programmingLanguage"] = "Nextflow DSL2"
                entity["license"] = {"@id": "https://www.apache.org/licenses/LICENSE-2.0"}
            graph.append(entity)
            workflow_dataset["hasPart"].append({"@id": identifier})

        params_id = "workflow/params.json"
        add_file(archive, args.params, params_id)
        graph.append(file_entity(params_id, args.params))
        workflow_dataset["hasPart"].append({"@id": params_id})

        archive.writestr(
            "ro-crate-metadata.json",
            json.dumps(metadata, indent=2, sort_keys=True) + "\n",
            compress_type=ZIP_DEFLATED,
        )

    print(output)


if __name__ == "__main__":
    main()
