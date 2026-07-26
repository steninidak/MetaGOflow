nextflow.enable.dsl = 2

def enabled(value) {
    if (value instanceof Boolean)
        return value
    return value?.toString()?.toBoolean()
}

include { QC }                    from './nextflow/workflows/qc'
include { MOTUS }                 from './nextflow/workflows/taxonomy'
include { RNA_TAXONOMY }          from './nextflow/workflows/taxonomy'
include { COMBINED_GENE_CALLER }  from './nextflow/workflows/gene_calling'
include { FUNCTIONAL_ANNOTATION } from './nextflow/workflows/functional_annotation'
include { MEGAHIT }               from './nextflow/workflows/assembly'

/*
 * MetaGOflow DSL2 entry point.
 *
 * Unlike the CWL workflow, disabled stages do not require pseudo files. A stage
 * consumes either the output of its predecessor or the corresponding
 * --processed_* restart parameter.
 */
workflow {
    if (params.help) {
        log.info """
        MetaGOflow (Nextflow)

        Required for a complete run:
          --reads 'sample_{1,2}.fastq.gz'
          --db_dir /path/to/ref-dbs

        Stage switches:
          --qc_and_merge_step          ${params.qc_and_merge_step}
          --taxonomic_inventory        ${params.taxonomic_inventory}
          --cgc_step                    ${params.cgc_step}
          --reads_functional_annotation ${params.reads_functional_annotation}
          --assemble                   ${params.assemble}

        Example:
          nextflow run main.nf -profile singularity \\
            --reads 'test_input/wgs-paired-SRR1620013_{1,2}.fastq.gz' \\
            --db_dir ref-dbs --outdir results
        """.stripIndent()
        return
    }

    Channel.empty().set { merged_ch }
    Channel.empty().set { paired_fasta_ch }
    Channel.empty().set { motus_input_ch }
    Channel.empty().set { ncrna_ch }
    Channel.empty().set { proteins_ch }

    if (enabled(params.qc_and_merge_step)) {
        if (!params.reads)
            error "QC is enabled: provide paired reads with --reads"
        raw_reads = Channel.fromFilePairs(params.reads, checkIfExists: true, flat: true)
        QC(raw_reads)
        merged_ch       = QC.out.merged_fasta
        paired_fasta_ch = QC.out.filtered_pair.collect()
        motus_input_ch  = QC.out.motus_input
    }
    else {
        if (params.processed_reads)
            merged_ch = Channel.value(file(params.processed_reads, checkIfExists: true))
        if (params.processed_read_files)
            paired_fasta_ch = Channel.value(params.processed_read_files.collect { file(it, checkIfExists: true) })
        if (params.input_for_motus)
            motus_input_ch = Channel.value(file(params.input_for_motus, checkIfExists: true))
    }

    if (enabled(params.taxonomic_inventory)) {
        if (!params.db_dir)
            error "Taxonomic inventory is enabled: provide --db_dir"
        motus_reads_ch = enabled(params.qc_and_merge_step)
            ? motus_input_ch
            : (params.input_for_motus ? motus_input_ch : merged_ch)
        MOTUS(motus_reads_ch)
        RNA_TAXONOMY(merged_ch)
        ncrna_ch = RNA_TAXONOMY.out.ncrna
    }
    else if (params.maskfile) {
        ncrna_ch = Channel.value(file(params.maskfile, checkIfExists: true))
    }

    if (enabled(params.cgc_step)) {
        COMBINED_GENE_CALLER(merged_ch)
        proteins_ch = COMBINED_GENE_CALLER.out.faa
    }
    else if (params.predicted_faa_from_previous_run) {
        proteins_ch = Channel.value(file(params.predicted_faa_from_previous_run, checkIfExists: true))
    }

    if (enabled(params.reads_functional_annotation)) {
        FUNCTIONAL_ANNOTATION(proteins_ch, ncrna_ch)
    }

    if (enabled(params.assemble)) {
        MEGAHIT(paired_fasta_ch)
    }
}
