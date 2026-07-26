process FASTP {
    tag "${sample_id}"
    label 'medium'
    container 'microbiomeinformatics/pipeline-v5.fastp:0.20.0'
    publishDir "${params.outdir}/qc", mode: params.publish_mode

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_1.trimmed.fastq"), path("${sample_id}_2.trimmed.fastq"), emit: trimmed
    tuple val(sample_id), path("${sample_id}.merged.fastq"), emit: merged
    path "${sample_id}.fastp.html", emit: html
    path "${sample_id}.fastp.json", emit: json

    script:
    def opts = [
        params.detect_adapter_for_pe ? '--detect_adapter_for_pe' : '',
        params.overrepresentation_analysis ? '--overrepresentation_analysis' : '',
        params.force_polyg_tail_trimming ? '--trim_poly_g' : '',
        params.disable_trim_poly_g ? '--disable_trim_poly_g' : '',
        params.base_correction ? '--correction' : '',
        params.cut_right ? '--cut_right' : '',
        params.qualified_phred_quality != null ? "--qualified_quality_phred ${params.qualified_phred_quality}" : '',
        params.unqualified_percent_limit != null ? "--unqualified_percent_limit ${params.unqualified_percent_limit}" : ''
    ].findAll().join(' ')
    """
    fastp \
      --in1 ${read1} --in2 ${read2} \
      --out1 ${sample_id}_1.trimmed.fastq \
      --out2 ${sample_id}_2.trimmed.fastq \
      --merge --merged_out ${sample_id}.merged.fastq \
      --length_required ${params.min_length_required} \
      --overlap_len_require ${params.overlap_len_require} \
      --thread ${task.cpus} \
      --html ${sample_id}.fastp.html --json ${sample_id}.fastp.json \
      ${opts}
    """
}

process FASTQ_TO_FILTERED_FASTA {
    tag "${sample_id}:${label}"
    label 'medium'
    container 'microbiomeinformatics/pipeline-v5.python2:v1'

    input:
    tuple val(sample_id), val(label), path(fastq)

    output:
    tuple val(sample_id), val(label), path("${sample_id}.${label}.fasta"), emit: fasta
    tuple val(sample_id), val(label), path("${sample_id}.${label}.qc_summary"), emit: summary

    script:
    """
    tr '\" /|<_;#' '-------' < ${fastq} > cleaned.fastq
    awk 'NR%4==1 {sub(/^@/,">"); print; next} NR%4==2 {print}' cleaned.fastq > unfiltered.fasta
    run_quality_filtering.py \
      unfiltered.fasta ${sample_id}.${label}.fasta ${sample_id}.${label}.qc_summary \
      \$(awk 'END {print int(NR/4)}' cleaned.fastq) \
      --min_length ${params.min_length_required} --extension fasta
    """
}

process QC_STATS {
    tag "${sample_id}:${label}"
    label 'light'
    container 'microbiomeinformatics/pipeline-v5.python2:v1'
    publishDir "${params.outdir}/qc/stats", mode: params.publish_mode

    input:
    tuple val(sample_id), val(label), path(fasta)

    output:
    tuple val(sample_id), val(label), path("${sample_id}.${label}.qc"), emit: stats

    script:
    """
    mkdir -p ${sample_id}.${label}.qc
    count=\$(grep -c '^>' ${fasta} || true)
    MGRAST_base.py \
      -o ${sample_id}.${label}.qc/summary \
      -d ${sample_id}.${label}.qc/nucleotide_distribution.full \
      -g ${sample_id}.${label}.qc/gc_sum.full \
      -l ${sample_id}.${label}.qc/length_sum.full \
      ${fasta}
    printf '%s\\n' "\$count" > ${sample_id}.${label}.qc/sequence_count
    touch ${sample_id}.${label}.qc/\$([ "\$count" -gt 0 ] && echo QC-PASSED || echo QC-FAILED)
    """
}

process INPUT_CHECKSUM {
    tag "${sample_id}"
    label 'light'
    container 'microbiomeinformatics/pipeline-v5.python3:v3.1'
    publishDir "${params.outdir}/qc/checksums", mode: params.publish_mode

    input:
    tuple val(sample_id), path(read)

    output:
    path "${read.simpleName}.sha1"

    script:
    """
    sha1sum ${read} > ${read.simpleName}.sha1
    """
}

workflow QC {
    take:
    reads

    main:
    FASTP(reads)

    FASTP.out.trimmed
        .flatMap { id, r1, r2 -> [[id, '1', r1], [id, '2', r2]] }
        .set { trimmed_reads }

    FASTP.out.merged
        .map { id, fq -> [id, 'merged', fq] }
        .mix(trimmed_reads)
        .set { all_fastq }

    FASTQ_TO_FILTERED_FASTA(all_fastq)
    QC_STATS(FASTQ_TO_FILTERED_FASTA.out.fasta)

    reads
        .flatMap { id, r1, r2 -> [[id, r1], [id, r2]] }
        .set { checksum_reads }
    INPUT_CHECKSUM(checksum_reads)

    emit:
    merged_fasta = FASTQ_TO_FILTERED_FASTA.out.fasta
        .filter { id, label, f -> label == 'merged' }
        .map { id, label, f -> f }
    motus_input = FASTP.out.merged.map { id, f -> f }
    filtered_pair = FASTQ_TO_FILTERED_FASTA.out.fasta
        .filter { id, label, f -> label in ['1','2'] }
        .map { id, label, f -> f }
    summaries = FASTQ_TO_FILTERED_FASTA.out.summary
    stats = QC_STATS.out.stats
}

