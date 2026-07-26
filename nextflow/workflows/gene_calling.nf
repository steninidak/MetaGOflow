process SPLIT_NUCLEOTIDES {
    tag "${fasta.simpleName}"
    label 'light'
    container 'microbiomeinformatics/pipeline-v5.split-fasta:v2'

    input:
    path fasta

    output:
    path 'chunks/*.fasta', emit: chunks

    script:
    """
    mkdir chunks
    split_fasta_by_size.sh -n ${params.cgc_chunk_size} ${fasta}
    mv *.fasta.* chunks/ 2>/dev/null || mv *_chunk* chunks/
    for f in chunks/*; do mv "\$f" "\${f}.fasta"; done
    """
}

process FRAGGENESCAN {
    tag "${chunk.simpleName}"
    label 'heavy'
    container 'hariszaf/pipeline-v5.fraggenescan:v1.31.1'

    input:
    path chunk

    output:
    tuple val(chunk.simpleName), path("${chunk.simpleName}.fgs.out"),
        path("${chunk.simpleName}.fgs.ffn"), path("${chunk.simpleName}.fgs.faa")

    script:
    """
    run_FGS.sh -i ${chunk} -o ${chunk.simpleName}.fgs -s 0 -t illumina_10
    """
}

process POSTPROCESS_GENES {
    tag "${name}"
    label 'medium'
    container 'microbiomeinformatics/pipeline-v5.protein-post-processing:v1.0.1'

    input:
    tuple val(name), path(fgs_out), path(fgs_ffn), path(fgs_faa)
    output:
    tuple val(name), path("${name}.faa"), path("${name}.ffn")

    script:
    """
    unite_protein_predictions.py \
      --fgs-out ${fgs_out} --fgs-ffn ${fgs_ffn} --fgs-faa ${fgs_faa} \
      --name ${name}
    """
}

process CONCATENATE_GENES {
    label 'light'
    publishDir "${params.outdir}/gene-calling", mode: params.publish_mode

    input:
    path faa_files
    path ffn_files

    output:
    path 'predicted_CDS.faa', emit: faa
    path 'predicted_CDS.ffn', emit: ffn
    path 'predicted_CDS.count', emit: count

    script:
    """
    cat ${faa_files} > predicted_CDS.faa
    cat ${ffn_files} > predicted_CDS.ffn
    grep -c '^>' predicted_CDS.faa > predicted_CDS.count || true
    """
}

workflow COMBINED_GENE_CALLER {
    take:
    fasta

    main:
    SPLIT_NUCLEOTIDES(fasta)
    FRAGGENESCAN(SPLIT_NUCLEOTIDES.out.chunks.flatten())
    POSTPROCESS_GENES(FRAGGENESCAN.out)
    POSTPROCESS_GENES.out
        .map { name, faa, ffn -> faa }
        .collect()
        .set { all_faa }
    POSTPROCESS_GENES.out
        .map { name, faa, ffn -> ffn }
        .collect()
        .set { all_ffn }
    CONCATENATE_GENES(all_faa, all_ffn)

    emit:
    faa = CONCATENATE_GENES.out.faa
    ffn = CONCATENATE_GENES.out.ffn
    count = CONCATENATE_GENES.out.count
}
