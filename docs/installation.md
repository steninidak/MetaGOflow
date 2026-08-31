# MetaGOflow Nextflow installation manual

This manual installs the Nextflow implementation of MetaGOflow on a Linux
workstation or a Slurm HPC system. Commands are run from the repository root,
the directory containing `main.nf`.

## 1. What is installed

MetaGOflow uses:

- Nextflow DSL2 to schedule tasks;
- Java to run Nextflow;
- Docker, Apptainer, or Singularity to provide analysis software;
- the MetaGOflow reference bundle for SILVA, Rfam, KOfam, eggNOG, and
  InterProScan;
- Python 3 from the host only for final RO-Crate creation.

Analysis programs such as fastp, mOTUs, Infernal, FragGeneScan, eggNOG-mapper,
InterProScan, HMMER, and MEGAHIT do not need to be installed on the host. They
are supplied by process-specific containers.

## 2. Hardware and storage planning

Use a Linux filesystem visible to every compute node that will execute jobs.
For a complete run, plan conservatively for:

- at least 300 GB for reference data;
- additional space for the container cache;
- a work directory several times larger than the compressed input;
- a second copy of published results when `publish_mode = 'copy'`;
- another result-sized allocation while the RO-Crate ZIP is being written.

Functional annotation, especially InterProScan and eggNOG, is the most
resource-intensive part. A small QC-only test can run on a workstation.
Production functional annotation is better suited to a scheduler-backed HPC
system.

Do not put `work/`, the reference databases, or the container cache on a small
home filesystem. On HPC, use project or scratch storage with sufficient inode
capacity.

## 3. Obtain the repository

```bash
git clone https://github.com/steninidak/MetaGOflow-Nextflow.git
cd MetaGOflow-Nextflow
```

The workflow entry point is `main.nf`. Process and subworkflow definitions are
kept under `nextflow/`, scheduler configuration templates under `conf/`, and
installation helpers under `Installation/`.

## 4. Install Java

Nextflow currently requires Bash 3.2 or newer and Java 17 or newer. Java 17 or
21 LTS is recommended.

Check an existing installation:

```bash
bash --version
java -version
```

On systems with environment modules:

```bash
module avail java
module load java/17
java -version
```

On Debian or Ubuntu, a system administrator can install OpenJDK:

```bash
sudo apt update
sudo apt install openjdk-17-jre-headless
```

If multiple Java installations exist, set `JAVA_HOME` to the selected Java
installation before running Nextflow. Consult the
[official Nextflow installation documentation](https://docs.seqera.io/nextflow/install)
for SDKMAN and other supported installation methods.

## 5. Install Nextflow

The recommended self-install method is:

```bash
curl -s https://get.nextflow.io | bash
chmod +x nextflow
mkdir -p "$HOME/.local/bin"
mv nextflow "$HOME/.local/bin/"
export PATH="$HOME/.local/bin:$PATH"
nextflow info
```

Add the `PATH` export to `.bashrc` or the site module configuration to make it
persistent. MetaGOflow requires Nextflow 24.04 or newer.

For an offline cluster, download a `nextflow-*-dist` standalone executable and
required container images on an internet-connected system, then transfer them
according to local policy. Note that Nextflow plugins may also require network
access unless already cached.

## 6. Choose one container runtime

Do not enable multiple container profiles in the same run.

### 6.1 Apptainer for HPC

Apptainer is recommended on shared clusters because it normally runs
containers without a privileged daemon. Ask the HPC administrators whether a
site module is already available:

```bash
module avail apptainer
module load apptainer
apptainer --version
apptainer exec docker://alpine:3.20 echo "Apptainer works"
```

If administrators need to install it, use the
[official Apptainer installation guide](https://apptainer.org/docs/admin/latest/installation.html).

MetaGOflow stores downloaded images in `.apptainer/` by default. Override this
with a site configuration if the repository filesystem is unsuitable:

```groovy
apptainer.cacheDir = '/project/my-group/container-cache'
```

### 6.2 Singularity

For sites that provide SingularityCE:

```bash
module load singularity
singularity --version
singularity exec docker://alpine:3.20 echo "Singularity works"
```

The pipeline uses `.singularity/` as its default cache directory.

### 6.3 Docker for a workstation

Install Docker Engine using the instructions for the host distribution in the
[official Docker installation guide](https://docs.docker.com/engine/install/).
After installation:

```bash
docker version
docker run --rm hello-world
```

If `docker version` reports a permission error, follow Docker's official
post-installation or rootless-mode guidance. Membership of the `docker` group
is effectively root-level access; use it only according to local security
policy.

## 7. Download the reference databases

The download is large and may take hours. Start it in a persistent terminal
such as `tmux`, or submit it as a transfer job if required by the site.

Install transfer utilities:

```bash
sudo apt install wget tar gzip
```

Download into a dedicated absolute path:

```bash
mkdir -p /project/my-group/metagoflow
bash Installation/download_dbs.sh \
  -f /project/my-group/metagoflow/ref-dbs
```

The downloader obtains fixed releases expected by this pipeline. Interrupted
downloads are not automatically transactional; inspect partial files before
restarting. Do not mix database releases silently, particularly InterProScan
software and data releases.

The expected layout is:

```text
ref-dbs/
├── Rfam/
│   ├── other/*.cm
│   ├── ribosomal/*.cm
│   └── rRNA.claninfo
├── db_kofam/
│   ├── db_kofam.hmm
│   └── db_kofam.hmm.h3{f,i,m,p}
├── eggnog/
│   ├── eggnog.db
│   └── eggnog_proteins.dmnd
├── interproscan-5.77-108.0/data/
├── kofam_ko_desc.tsv
├── silva_lsu/
│   ├── LSU.fasta
│   ├── LSU.fasta.mscluster
│   ├── lsu2.otu
│   └── slv_lsu_filtered2.txt
└── silva_ssu/
    ├── SSU.fasta
    ├── SSU.fasta.mscluster
    ├── ssu2.otu
    └── slv_ssu_filtered2.txt
```

If reference data is managed centrally or has a different release name,
override functional database locations:

```bash
--interproscan_data_dir /refs/interproscan-5.77-108.0/data \
--eggnog_data_dir /refs/eggnog \
--hmm_database /refs/db_kofam/db_kofam.hmm \
--hmm_database_dir /refs/db_kofam \
--ko_description_file /refs/kofam_ko_desc.tsv
```

## 8. Validate the installation

From the repository root, make the checker executable if needed:

```bash
chmod +x Installation/check_nextflow_install.sh
```

Check taxonomy dependencies with Apptainer:

```bash
Installation/check_nextflow_install.sh \
  -d /project/my-group/metagoflow/ref-dbs \
  -p apptainer
```

Include functional-annotation databases:

```bash
Installation/check_nextflow_install.sh \
  -d /project/my-group/metagoflow/ref-dbs \
  -p apptainer -f
```

For Docker, use `-p docker`; for Singularity, use `-p singularity`. Fix every
`ERROR` before a production run.

## 9. First execution

### 9.1 QC-only smoke test

This confirms Nextflow, the selected runtime, input pairing, the fastp image,
publishing, reporting, and RO-Crate packaging without reference databases:

```bash
nextflow run main.nf -profile apptainer \
  --reads 'test_input/wgs-paired-SRR1620013_{1,2}.fastq.gz' \
  --taxonomic_inventory false \
  --cgc_step false \
  --reads_functional_annotation false \
  --assemble false \
  --outdir results-smoke
```

Expected top-level products include:

```text
results-smoke/
├── pipeline_info/
├── qc/
└── ro-crate/<run-name>.ro-crate.zip
```

Inspect the execution trace and verify the crate:

```bash
column -t results-smoke/pipeline_info/trace.txt | less -S
unzip -t results-smoke/ro-crate/*.ro-crate.zip
```

### 9.2 Taxonomy run

```bash
nextflow run main.nf -profile apptainer \
  --reads 'test_input/wgs-paired-SRR1620013_{1,2}.fastq.gz' \
  --db_dir /project/my-group/metagoflow/ref-dbs \
  --qc_and_merge_step true \
  --taxonomic_inventory true \
  --cgc_step false \
  --reads_functional_annotation false \
  --assemble false \
  --outdir results-taxonomy
```

### 9.3 Complete run

```bash
nextflow run main.nf -profile apptainer \
  --reads '/data/sample_{1,2}.fastq.gz' \
  --db_dir /project/my-group/metagoflow/ref-dbs \
  --qc_and_merge_step true \
  --taxonomic_inventory true \
  --cgc_step true \
  --reads_functional_annotation true \
  --assemble true \
  --threads 32 \
  --interproscan_threads 8 \
  --max_cpus 64 \
  --max_memory '256 GB' \
  --outdir results/sample
```

Use quotes around brace-based read patterns so the shell does not expand them
before Nextflow. The pattern must resolve to exactly two mates per sample.

## 10. Slurm installation and configuration

Nextflow must run on a login or submission node where `sbatch`, `squeue`, and
`scancel` are available. The launch directory, work directory, reference
directory, input files, and container cache must be visible from compute
nodes.

Copy the template instead of editing the tracked file:

```bash
cp conf/slurm.config conf/site-slurm.config
```

Example:

```groovy
process {
    executor = 'slurm'
    queue = 'compute'
    clusterOptions = '--account=my-project'
}

executor {
    queueSize = 100
    submitRateLimit = '20 sec'
}

apptainer.cacheDir = '/project/my-project/apptainer-cache'
```

Launch:

```bash
nextflow run main.nf \
  -profile apptainer,slurm \
  -c conf/site-slurm.config \
  -work-dir /scratch/my-project/metagoflow-work \
  --reads '/data/sample_{1,2}.fastq.gz' \
  --db_dir /project/my-project/ref-dbs \
  --outdir /project/my-project/results/sample
```

Nextflow submits each process as a separate Slurm job. Do not wrap the entire
pipeline in a large compute-node allocation unless required by site policy.
See the [official Nextflow executor documentation](https://docs.seqera.io/nextflow/executor)
for scheduler semantics.

## 11. Resume, retry, and upgrades

Resume an interrupted run by repeating the same command with:

```bash
-resume
```

Keep the `work/` directory and `.nextflow/` cache until results have been
accepted. Changing input content, commands, relevant parameters, containers,
or process configuration invalidates affected cached tasks.

Before upgrading Nextflow or changing container/database releases:

1. preserve the command, parameter file, trace, and RO-Crate;
2. test the QC-only sample;
3. run a representative taxonomy or functional sample;
4. compare counts and expected output files;
5. only then upgrade production runs.

## 12. Common installation failures

### Java is too old

Symptom: Nextflow reports that Java is missing or unsupported.

```bash
java -version
echo "$JAVA_HOME"
```

Load or install Java 17 or 21 and ensure `java` resolves to it.

### Container command not found

The selected profile must match the installed runtime:

```bash
-profile docker
-profile apptainer
-profile singularity
```

### Container downloads fail on compute nodes

Pre-populate a shared cache from a node with registry access, or ask the
administrator to mirror the images. Confirm that compute nodes can read the
cache.

### Docker permission denied

Confirm `docker run --rm hello-world` works as the same user launching
Nextflow. Do not work around daemon permissions with arbitrary `sudo nextflow`
commands because this creates root-owned work files.

### Missing InterProScan data

Set `--interproscan_data_dir` to the actual `data/` directory. The data release
must be compatible with the InterProScan container.

### No space left on device

Check both bytes and inodes:

```bash
df -h .
df -ih .
du -sh work results* .apptainer .singularity 2>/dev/null
```

Move `-work-dir` and the container cache to suitable storage. Use
`nextflow clean` only after confirming that resumption is no longer required.

### RO-Crate fails or is unexpectedly large

RO-Crate packaging includes result payloads and calculates SHA-256 checksums.
It therefore needs time and temporary space proportional to the results.
Disable it for a diagnostic run with `--ro_crate false`; leave it enabled for
archival production runs.

## 13. Installation acceptance checklist

- `java -version` reports Java 17 or newer.
- `nextflow info` succeeds.
- Exactly one chosen container runtime executes a test image.
- The database checker completes without errors for enabled stages.
- Reference and work paths are visible on compute nodes.
- The QC-only smoke test succeeds.
- `pipeline_info/trace.txt` contains successful processes.
- The RO-Crate ZIP passes `unzip -t`.
- A representative biological result has been reviewed before production use.
