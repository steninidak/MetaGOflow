process RUN_MOTUS {
    tag "${reads.simpleName}"
    label 'heavy'
    container 'microbiomeinformatics/pipeline-v5.motus:v2.5.1'
    publishDir "${params.outdir}/taxonomy/motus", mode: params.publish_mode

    input:
    path reads

    output:
    path "${reads.simpleName}.motus.tsv", emit: taxonomy

    script:
    """
    motus profile -c -q -s ${reads} -t ${task.cpus} > raw.motus.tsv
    clean_motus_output.sh raw.motus.tsv > ${reads.simpleName}.motus.tsv
    """
}

process CMSEARCH {
    tag "${model.simpleName}"
    label 'heavy'
    container 'microbiomeinformatics/pipeline-v5.cmsearch:v1.1.2'

    input:
    tuple path(reads), path(model)

    output:
    tuple val(model.simpleName), path("${reads.simpleName}.${model.simpleName}.cmsearch.tbl")

    script:
    """
    cmsearch --cpu ${task.cpus} --cut_ga --noali \
      --tblout ${reads.simpleName}.${model.simpleName}.cmsearch.tbl \
      ${model} ${reads} > ${reads.simpleName}.${model.simpleName}.cmsearch.log
    """
}

process DEOVERLAP_RNA {
    tag "${reads.simpleName}"
    label 'medium'
    container 'microbiomeinformatics/pipeline-v5.cmsearch-deoverlap:v0.02'

    input:
    path reads
    path tables
    path clan_info

    output:
    path "${reads.simpleName}.cmsearch.all.tblout", emit: all_hits
    path "${reads.simpleName}.cmsearch.all.tblout.deoverlapped", emit: ncrna

    script:
    """
    cat ${tables} > ${reads.simpleName}.cmsearch.all.tblout
    cmsearch-deoverlap.pl \
      --clanin ${clan_info} \
      --tblout ${reads.simpleName}.cmsearch.all.tblout.deoverlapped \
      ${reads.simpleName}.cmsearch.all.tblout
    """
}

process EXTRACT_RNA_SEQUENCES {
    tag "${reads.simpleName}"
    label 'medium'
    container 'microbiomeinformatics/pipeline-v5.easel:v0.45h'

    input:
    path reads
    path ncrna

    output:
    tuple val(reads.simpleName), path('extracted.fasta')

    script:
    """
    awk '{print \$1"-"\$3"/q"\$8"-"\$9" "\$8" "\$9" "\$1}' ${ncrna} > coords.txt
    cp ${reads} indexed.fasta
    esl-sfetch --index indexed.fasta
    esl-sfetch -Cf indexed.fasta coords.txt > extracted.fasta
    """
}

process CATEGORIZE_RNA {
    tag "${sample_id}"
    label 'medium'
    container 'microbiomeinformatics/pipeline-v5.python3:v3.1'

    input:
    tuple val(sample_id), path(extracted)

    output:
    path "${sample_id}.SSU.fasta", emit: ssu
    path "${sample_id}.LSU.fasta", emit: lsu
    path "${sample_id}.rna-count.tsv", emit: counts
    path "${sample_id}.*.fasta.gz", emit: categories

    script:
    """
    get_subunits.py -i ${extracted} \
      -s '${params.ssu_label}' -l '${params.lsu_label}' \
      -f '${params.pattern_5s}' -e '${params.pattern_5_8s}' -p '${sample_id}'
    find sequence-categorisation -type f -name '*SSU.fasta*' -exec cp {} ${sample_id}.SSU.fasta \\;
    find sequence-categorisation -type f -name '*LSU.fasta*' -exec cp {} ${sample_id}.LSU.fasta \\;
    test -f ${sample_id}.SSU.fasta || touch ${sample_id}.SSU.fasta
    test -f ${sample_id}.LSU.fasta || touch ${sample_id}.LSU.fasta
    printf 'SSU\\t%s\\nLSU\\t%s\\n' \
      "\$(grep -c '^>' ${sample_id}.SSU.fasta || true)" \
      "\$(grep -c '^>' ${sample_id}.LSU.fasta || true)" \
      > ${sample_id}.rna-count.tsv
    cp sequence-categorisation/*.fa . 2>/dev/null || true
    for fasta in *.fa ${sample_id}.*.fasta; do
      test -f "\$fasta" || continue
      gzip -c "\$fasta" > "${sample_id}.\$(basename "\$fasta").gz"
    done
    """
}

process CLASSIFY_RRNA {
    tag "${label}"
    label 'heavy'
    container 'microbiomeinformatics/pipeline-v5.mapseq:v1.2.3'
    publishDir { "${params.outdir}/taxonomy/${label}" }, mode: params.publish_mode

    input:
    tuple val(label), path(fasta), path(database), path(taxonomy), path(otus)

    output:
    path "${label}", emit: directory

    script:
    """
    mkdir -p ${label}
    count=\$(grep -c '^>' ${fasta} || true)
    if [ "\$count" -gt 0 ]; then
      mapseq -nthreads ${task.cpus} -tophits 80 -topotus 40 -outfmt simple \
        -db ${database} -tax ${taxonomy} ${fasta} > ${label}/mapseq.tsv
      mapseq2biom.pl --otu ${otus} --query ${fasta} \
        --label ${label} --outfile ${label}/otu.tsv \
        --krona ${label}/krona.txt --notaxidfile ${label}/otu.notaxid.tsv \
        < ${label}/mapseq.tsv
      ktImportText -o ${label}/krona.html ${label}/krona.txt
    else
      touch ${label}/mapseq.tsv ${label}/otu.tsv ${label}/krona.txt ${label}/krona.html
    fi
    """
}

process PACKAGE_TAXONOMY {
    label 'light'
    publishDir "${params.outdir}/taxonomy", mode: params.publish_mode

    input:
    path category_files
    path tax_dirs

    output:
    path 'sequence-categorisation', emit: sequence_categories
    path 'taxonomy-summary', emit: taxonomy_summary

    script:
    """
    mkdir sequence-categorisation taxonomy-summary
    cp -L ${category_files} sequence-categorisation/
    cp -RL ${tax_dirs} taxonomy-summary/
    """
}

workflow MOTUS {
    take:
    reads
    main:
    RUN_MOTUS(reads)
    emit:
    taxonomy = RUN_MOTUS.out.taxonomy
}

workflow RNA_TAXONOMY {
    take:
    reads

    main:
    db_root = file(params.db_dir, checkIfExists: true)
    models = Channel.fromPath("${params.db_dir}/Rfam/**/*.cm", checkIfExists: true)
    reads.combine(models).map { r, m -> tuple(r, m) }.set { searches }
    CMSEARCH(searches)

    clan = Channel.value(file("${params.db_dir}/Rfam/rRNA.claninfo", checkIfExists: true))
    hit_tables = CMSEARCH.out.collect { name, table -> table }
    DEOVERLAP_RNA(reads, hit_tables, clan)
    EXTRACT_RNA_SEQUENCES(reads, DEOVERLAP_RNA.out.ncrna)
    CATEGORIZE_RNA(EXTRACT_RNA_SEQUENCES.out)

    ssu_refs = Channel.value(tuple(
        'SSU',
        file("${params.db_dir}/silva_ssu/SSU.fasta", checkIfExists: true),
        file("${params.db_dir}/silva_ssu/slv_ssu_filtered2.txt", checkIfExists: true),
        file("${params.db_dir}/silva_ssu/ssu2.otu", checkIfExists: true)
    ))
    lsu_refs = Channel.value(tuple(
        'LSU',
        file("${params.db_dir}/silva_lsu/LSU.fasta", checkIfExists: true),
        file("${params.db_dir}/silva_lsu/slv_lsu_filtered2.txt", checkIfExists: true),
        file("${params.db_dir}/silva_lsu/lsu2.otu", checkIfExists: true)
    ))
    CATEGORIZE_RNA.out.ssu.combine(ssu_refs).map { f, l, d, t, o -> tuple(l,f,d,t,o) }
        .mix(CATEGORIZE_RNA.out.lsu.combine(lsu_refs).map { f,l,d,t,o -> tuple(l,f,d,t,o) })
        .set { classify_inputs }
    CLASSIFY_RRNA(classify_inputs)
    PACKAGE_TAXONOMY(CATEGORIZE_RNA.out.categories.collect(), CLASSIFY_RRNA.out.directory.collect())

    emit:
    ncrna = DEOVERLAP_RNA.out.ncrna
    counts = CATEGORIZE_RNA.out.counts
    sequence_categories = PACKAGE_TAXONOMY.out.sequence_categories
    taxonomy_summary = PACKAGE_TAXONOMY.out.taxonomy_summary
}
