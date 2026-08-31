# Running MetaGOflow with Nextflow

MetaGOflow is implemented in Nextflow DSL2 with containerized tools, reference
databases, five stage switches, and intermediate-file entry points. Nextflow
provides resumable content-addressed tasks, native scatter through channels,
and executor-specific scheduling.

## Requirements

- Nextflow 24.04 or newer
- Java 17 or newer
- Docker, Singularity, or Apptainer
- The existing MetaGOflow reference database bundle

## Complete run

```bash
nextflow run main.nf -profile singularity \
  --reads 'test_input/wgs-paired-SRR1620013_{1,2}.fastq.gz' \
  --db_dir ref-dbs \
  --qc_and_merge_step true \
  --taxonomic_inventory true \
  --cgc_step true \
  --reads_functional_annotation true \
  --assemble true \
  --outdir results
```

Use `-profile docker` or `-profile apptainer` for the other supported
container runtimes. Add `-profile singularity,slurm` on a Slurm cluster. Site
settings such as account, partition, and cluster options belong in a local
configuration file:

```bash
nextflow run main.nf -profile singularity,slurm \
  -c conf/slurm.config -resume [parameters...]
```

## Resume and partial runs

Retry an interrupted run by repeating the command with `-resume`. Nextflow
reuses successful tasks whose inputs, commands, containers, and relevant
configuration have not changed.

Disabled upstream stages accept real intermediate files:

```bash
nextflow run main.nf -profile singularity \
  --qc_and_merge_step false \
  --taxonomic_inventory false \
  --cgc_step false \
  --reads_functional_annotation true \
  --processed_reads results/sample.merged.fasta \
  --maskfile results/sample.cmsearch.all.tblout.deoverlapped \
  --predicted_faa_from_previous_run results/sample_CDS.faa \
  --db_dir ref-dbs --outdir functional-rerun
```

For assembly without QC, pass a JSON/Groovy list:

```bash
--processed_read_files '["forward.fasta","reverse.fasta"]'
```

## Resource parameters

Global resource ceilings are `max_cpus`, `max_memory`, and `max_time`. Per-site
overrides can be assigned by process label (`light`, `medium`, or `heavy`) or
by process name in a Nextflow configuration file.

## Outputs and provenance

Outputs are grouped beneath `outdir` by stage. `pipeline_info` contains the
execution trace, report, timeline, and DAG. The work directory is deliberately
kept outside the result set so `-resume` remains available; remove it only
after the run is accepted.

After a successful run, MetaGOflow creates
`outdir/ro-crate/<run-name>.ro-crate.zip`. The RO-Crate 1.1 archive contains
the analysis results, available execution trace, Nextflow sources, resolved
parameters, SHA-256 checksums, file media types, and workflow-run provenance.
The final HTML report and timeline remain in `outdir/pipeline_info`. Result
files remain in `outdir`, preserving Nextflow's `-resume` behavior. Disable
packaging with `--ro_crate false`.
