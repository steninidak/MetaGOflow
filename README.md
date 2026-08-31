# MetaGOflow

MetaGOflow is a Nextflow DSL2 workflow for marine Genomic Observatories'
metagenomic data analysis. It provides quality control, taxonomic inventory,
gene calling, functional annotation, assembly, and RO-Crate provenance.

## Requirements

- Nextflow 24.04 or newer
- Java 17 or newer
- Docker, Singularity, or Apptainer
- Approximately 300 GB for the reference database bundle

Analysis tools are supplied through process-specific containers and do not
need to be installed on the host.

## Installation

Clone the repository and download the databases:

```bash
git clone https://github.com/steninidak/MetaGOflow-Nextflow.git
cd MetaGOflow-Nextflow
bash Installation/download_dbs.sh -f ref-dbs
```

For complete host, container, database, and HPC setup instructions, see the
[installation manual](docs/installation.md).

## Run MetaGOflow

Run the complete workflow with Apptainer:

```bash
nextflow run main.nf -profile apptainer \
  --reads 'test_input/wgs-paired-SRR1620013_{1,2}.fastq.gz' \
  --db_dir ref-dbs \
  --qc_and_merge_step true \
  --taxonomic_inventory true \
  --cgc_step true \
  --reads_functional_annotation true \
  --assemble true \
  --outdir results
```

Use `-profile docker` or `-profile singularity` for another container runtime.
On Slurm, combine the runtime and scheduler profiles:

```bash
nextflow run main.nf -profile apptainer,slurm \
  -c conf/site-slurm.config \
  -work-dir /scratch/my-project/metagoflow-work \
  --reads '/data/sample_{1,2}.fastq.gz' \
  --db_dir /project/my-project/ref-dbs \
  --outdir /project/my-project/results/sample \
  -resume
```

See the [Nextflow guide](docs/nextflow.md) for complete and partial runs,
restart inputs, parameters, outputs, and cluster configuration.

## Test data

Small paired-read datasets are provided in `test_input/`. A QC-only smoke test
can be run with:

```bash
nextflow run main.nf -profile apptainer,test --db_dir ref-dbs
```

## Outputs

Results are grouped beneath `--outdir`. Execution reports, trace, timeline, and
DAG are written to `pipeline_info/`. A successful analysis also creates an
RO-Crate archive under `ro-crate/`. Keep the Nextflow work directory and
`.nextflow/` cache until the run has been accepted so that `-resume` remains
available.

## License

See [LICENSE](LICENSE).
