process RUN_MEGAHIT {
    tag 'assembly'
    label 'heavy'
    container 'quay.io/biocontainers/megahit:1.2.9--h2e03b76_1'
    publishDir "${params.outdir}/assembly", mode: params.publish_mode

    input:
    path read_pair

    output:
    path 'final.contigs.fa', emit: contigs
    path 'megahit.options.txt', emit: options

    script:
    def r1 = read_pair[0]
    def r2 = read_pair[1]
    """
    megahit -1 ${r1} -2 ${r2} \
      --min-contig-len ${params.min_contig_len} \
      --memory ${params.megahit_memory} \
      --num-cpu-threads ${task.cpus} \
      --out-dir megahit
    cp megahit/final.contigs.fa final.contigs.fa
    printf '%s\\n' \
      '--min-contig-len ${params.min_contig_len}' \
      '--memory ${params.megahit_memory}' \
      '--num-cpu-threads ${task.cpus}' > megahit.options.txt
    """
}

workflow MEGAHIT {
    take:
    reads
    main:
    RUN_MEGAHIT(reads)
    emit:
    contigs = RUN_MEGAHIT.out.contigs
    options = RUN_MEGAHIT.out.options
}

