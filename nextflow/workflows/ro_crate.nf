process CREATE_RO_CRATE {
    tag "${workflow.runName}"
    label 'light'
    executor 'local'
    publishDir "${params.outdir}/ro-crate", mode: params.publish_mode

    input:
    val completion_tokens
    path generator

    output:
    path "${workflow.runName}.ro-crate.zip", emit: archive

    script:
    def parameterJson = groovy.json.JsonOutput.prettyPrint(
        groovy.json.JsonOutput.toJson(params as Map)
    )
    def parameterBase64 = parameterJson.bytes.encodeBase64().toString()
    """
    printf '%s' '${parameterBase64}' | base64 --decode > params.json

    python3 ${generator} \
      --results '${file(params.outdir).toAbsolutePath()}' \
      --workflow '${projectDir}' \
      --params params.json \
      --output '${workflow.runName}.ro-crate.zip' \
      --run-name '${workflow.runName}' \
      --run-id '${workflow.sessionId}' \
      --nextflow-version '${workflow.nextflow.version}'
    """
}

workflow RO_CRATE {
    take:
    completion

    main:
    generator = channel.value(
        file("${projectDir}/nextflow/bin/create_ro_crate.py", checkIfExists: true)
    )
    CREATE_RO_CRATE(completion.collect(), generator)

    emit:
    archive = CREATE_RO_CRATE.out.archive
}
