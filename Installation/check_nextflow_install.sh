#!/usr/bin/env bash

set -uo pipefail

DB_DIR=""
PROFILE="apptainer"
CHECK_FUNCTIONAL=false

usage() {
    cat <<'EOF'
Usage: Installation/check_nextflow_install.sh -d DB_DIR [-p PROFILE] [-f]

Checks the host runtime, input examples, and reference files required by the
Nextflow pipeline. Use -f to include functional-annotation databases.

Options:
  -d DIR       Reference database directory
  -p PROFILE   docker, apptainer, or singularity (default: apptainer)
  -f           Check functional-annotation databases
  -h           Show this help
EOF
}

while getopts ":d:p:fh" option; do
    case "${option}" in
        d) DB_DIR=${OPTARG} ;;
        p) PROFILE=${OPTARG} ;;
        f) CHECK_FUNCTIONAL=true ;;
        h) usage; exit 0 ;;
        :) echo "Option -${OPTARG} requires a value" >&2; exit 2 ;;
        \?) echo "Unknown option -${OPTARG}" >&2; usage; exit 2 ;;
    esac
done

if [[ -z "${DB_DIR}" ]]; then
    echo "ERROR: -d DB_DIR is required" >&2
    usage
    exit 2
fi

failures=0

ok() { printf 'OK    %s\n' "$1"; }
bad() { printf 'ERROR %s\n' "$1" >&2; failures=$((failures + 1)); }
warn() { printf 'WARN  %s\n' "$1" >&2; }

check_command() {
    if command -v "$1" >/dev/null 2>&1; then
        ok "$1: $(command -v "$1")"
    else
        bad "command not found: $1"
    fi
}

check_file() {
    if [[ -s "$1" ]]; then
        ok "$1"
    else
        bad "missing or empty file: $1"
    fi
}

check_dir() {
    if [[ -d "$1" ]] && find "$1" -mindepth 1 -print -quit | grep -q .; then
        ok "$1"
    else
        bad "missing or empty directory: $1"
    fi
}

check_command java
check_command nextflow

java_major=$(java -version 2>&1 | awk -F '[".]' '/version/ {print ($2 == "1" ? $3 : $2); exit}')
if [[ "${java_major:-0}" =~ ^[0-9]+$ ]] && ((java_major >= 17)); then
    ok "Java major version ${java_major}"
else
    bad "Java 17 or newer is required; detected '${java_major:-unknown}'"
fi

case "${PROFILE}" in
    docker)
        check_command docker
        if docker info >/dev/null 2>&1; then
            ok "Docker daemon is reachable"
        else
            bad "Docker daemon is not reachable by this user"
        fi
        ;;
    apptainer)
        check_command apptainer
        ;;
    singularity)
        check_command singularity
        ;;
    *)
        bad "unsupported profile '${PROFILE}'"
        ;;
esac

check_file "main.nf"
check_file "nextflow.config"
check_file "test_input/wgs-paired-SRR1620013_1.fastq.gz"
check_file "test_input/wgs-paired-SRR1620013_2.fastq.gz"

check_file "${DB_DIR}/silva_ssu/SSU.fasta"
check_file "${DB_DIR}/silva_ssu/SSU.fasta.mscluster"
check_file "${DB_DIR}/silva_ssu/slv_ssu_filtered2.txt"
check_file "${DB_DIR}/silva_ssu/ssu2.otu"
check_file "${DB_DIR}/silva_lsu/LSU.fasta"
check_file "${DB_DIR}/silva_lsu/LSU.fasta.mscluster"
check_file "${DB_DIR}/silva_lsu/slv_lsu_filtered2.txt"
check_file "${DB_DIR}/silva_lsu/lsu2.otu"
check_file "${DB_DIR}/Rfam/rRNA.claninfo"
check_dir "${DB_DIR}/Rfam/ribosomal"
check_dir "${DB_DIR}/Rfam/other"

if [[ "${CHECK_FUNCTIONAL}" == true ]]; then
    check_file "${DB_DIR}/eggnog/eggnog.db"
    check_file "${DB_DIR}/eggnog/eggnog_proteins.dmnd"
    check_file "${DB_DIR}/db_kofam/db_kofam.hmm"
    check_file "${DB_DIR}/db_kofam/db_kofam.hmm.h3f"
    check_file "${DB_DIR}/db_kofam/db_kofam.hmm.h3i"
    check_file "${DB_DIR}/db_kofam/db_kofam.hmm.h3m"
    check_file "${DB_DIR}/db_kofam/db_kofam.hmm.h3p"
    check_file "${DB_DIR}/kofam_ko_desc.tsv"

    ips_count=$(find "${DB_DIR}" -maxdepth 2 -type d -path '*/interproscan-*/data' 2>/dev/null | wc -l)
    if ((ips_count > 0)); then
        ok "InterProScan data directory found"
    else
        bad "no DB_DIR/interproscan-*/data directory found"
    fi
fi

if ((failures > 0)); then
    printf '\nInstallation check failed with %d error(s).\n' "${failures}" >&2
    exit 1
fi

printf '\nInstallation check passed.\n'
