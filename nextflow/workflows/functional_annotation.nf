process SPLIT_PROTEINS {
    tag "${kind}:${proteins.simpleName}"
    label 'light'
    container 'microbiomeinformatics/pipeline-v5.python3:v3.1'

    input:
    tuple val(kind), val(size), path(proteins)

    output:
    tuple val(kind), path('chunks/*')

    script:
    """
    mkdir chunks
    cd chunks
    split_to_chunks.py -i ../${proteins} -s ${size} -f fasta
    """
}

process EGGNOG {
    tag "${chunk.simpleName}"
    label 'heavy'
    container 'cymon/eggnog-2.1.12:0.2'

    input:
    path chunk
    path eggnog_db
    path diamond_db
    path data_dir

    output:
    path "${chunk.simpleName}.emapper.annotations", emit: annotations
    path "${chunk.simpleName}.emapper.orthologs", optional: true, emit: orthologs

    script:
    """
    emapper.py -i ${chunk} -o ${chunk.simpleName} \
      --database ${eggnog_db} --dmnd_db ${diamond_db} \
      --data_dir ${data_dir} --cpu ${task.cpus}
    """
}

process INTERPROSCAN {
    tag "${chunk.simpleName}"
    label 'heavy'
    cpus { params.interproscan_threads }
    container 'sninidakis/ips:latest'

    input:
    path chunk
    path databases

    output:
    path "${chunk.simpleName}.IPS.tsv"

    script:
    def apps = params.interproscan_applications.join(',')
    """
    interproscan.sh --disable-precalc --goterms --pathways \
      --cpu ${task.cpus} --tempdir "\$PWD/tmp" \
      --input ${chunk} --applications '${apps}' -f TSV \
      --data-dir ${databases} -o ${chunk.simpleName}.IPS.tsv
    """
}

process HMMSEARCH {
    tag "${chunk.simpleName}"
    label 'heavy'
    container 'cymon/hmmer-3.4:0.1'

    input:
    path chunk
    path hmm_db
    path hmm_dir

    output:
    path "${chunk.simpleName}_hmmsearch.tbl"

    script:
    def noali = params.hmm_omit_alignment ? '--noali' : ''
    def cutga = params.hmm_gathering_bit_score ? '--cut_ga' : ''
    """
    hmmsearch ${noali} ${cutga} --cpu ${task.cpus} \
      --domtblout ${chunk.simpleName}_hmmsearch.tbl \
      -o /dev/null ${hmm_db} ${chunk}
    """
}

process COMBINE_FUNCTIONAL_RESULTS {
    label 'light'
    container 'microbiomeinformatics/pipeline-v5.python3:v3.1'
    publishDir "${params.outdir}/functional-annotation", mode: params.publish_mode

    input:
    path eggnog
    path ips
    path hmm

    output:
    path 'eggnog.annotations.tsv', emit: eggnog
    path 'interproscan.tsv', emit: ips
    path 'hmmsearch.tsv', emit: hmm

    script:
    """
    cat ${eggnog} > eggnog.annotations.tsv
    cat ${ips} > interproscan.tsv
    cat ${hmm} > hmmsearch.raw
    hmmscan_tab.py -i hmmsearch.raw -o hmmsearch.tsv
    """
}

process FUNCTIONAL_SUMMARIES {
    label 'medium'
    container 'microbiomeinformatics/pipeline-v5.python3:v3.1'
    publishDir "${params.outdir}/functional-annotation", mode: params.publish_mode

    input:
    path eggnog
    path ips
    path hmm
    path ncrna
    path proteins
    path ko_file

    output:
    path 'stats', emit: stats
    path 'summaries', emit: summaries

    script:
    """
    mkdir stats summaries
    awk -F '\\t' '\$4=="Pfam"' ${ips} > summaries/pfam.tsv
    awk -F '\\t' '!/^#/ && NF {print \$1"\\t"\$2}' ${eggnog} \
      | sort -u > summaries/eggnog.tsv
    functional_stats.py \
      --interproscan ${ips} --hmmscan ${hmm} \
      --pfam summaries/pfam.tsv --cmsearch_file ${ncrna} \
      --cds_file ${proteins} --ko_file ${ko_file} \
      > stats/functional_stats.txt
    printf 'artifact\\tlines\\n' > stats/counts.tsv
    for f in ${eggnog} ${ips} ${hmm} summaries/pfam.tsv; do
      printf '%s\\t%s\\n' "\$(basename "\$f")" "\$(wc -l < "\$f")" >> stats/counts.tsv
    done
    """
}

workflow FUNCTIONAL_ANNOTATION {
    take:
    proteins
    ncrna

    main:
    protein_splits = proteins
        .flatMap { p -> [
            tuple('eggnog', params.protein_chunk_size_eggnog as int, p),
            tuple('ips', params.protein_chunk_size_IPS as int, p),
            tuple('hmm', params.protein_chunk_size_hmm as int, p)
        ] }
    SPLIT_PROTEINS(protein_splits)

    egg_chunks = SPLIT_PROTEINS.out.filter { kind, f -> kind == 'eggnog' }.map { kind, f -> f }.flatten()
    ips_chunks = SPLIT_PROTEINS.out.filter { kind, f -> kind == 'ips' }.map { kind, f -> f }.flatten()
    hmm_chunks = SPLIT_PROTEINS.out.filter { kind, f -> kind == 'hmm' }.map { kind, f -> f }.flatten()

    egg_db = Channel.value(file("${params.db_dir}/eggnog/eggnog.db", checkIfExists: true))
    diamond_db = Channel.value(file("${params.db_dir}/eggnog/eggnog_proteins.dmnd", checkIfExists: true))
    egg_dir = Channel.value(file("${params.db_dir}/eggnog", checkIfExists: true))
    ips_dir = Channel.value(file("${params.db_dir}/interproscan-5.77-100.0/data", checkIfExists: true))
    hmm_db = Channel.value(file("${params.db_dir}/db_kofam/db_kofam.hmm", checkIfExists: true))
    hmm_dir = Channel.value(file("${params.db_dir}/db_kofam", checkIfExists: true))

    EGGNOG(egg_chunks, egg_db, diamond_db, egg_dir)
    INTERPROSCAN(ips_chunks, ips_dir)
    HMMSEARCH(hmm_chunks, hmm_db, hmm_dir)
    COMBINE_FUNCTIONAL_RESULTS(
        EGGNOG.out.annotations.collect(),
        INTERPROSCAN.out.collect(),
        HMMSEARCH.out.collect()
    )
    ko = Channel.value(file("${params.db_dir}/kofam_ko_desc.tsv", checkIfExists: true))
    FUNCTIONAL_SUMMARIES(
        COMBINE_FUNCTIONAL_RESULTS.out.eggnog,
        COMBINE_FUNCTIONAL_RESULTS.out.ips,
        COMBINE_FUNCTIONAL_RESULTS.out.hmm,
        ncrna, proteins, ko
    )

    emit:
    eggnog = COMBINE_FUNCTIONAL_RESULTS.out.eggnog
    ips = COMBINE_FUNCTIONAL_RESULTS.out.ips
    hmm = COMBINE_FUNCTIONAL_RESULTS.out.hmm
    stats = FUNCTIONAL_SUMMARIES.out.stats
    summaries = FUNCTIONAL_SUMMARIES.out.summaries
}

