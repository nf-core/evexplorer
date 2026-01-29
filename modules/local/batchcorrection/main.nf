process BATCHCORRECTION {
    // tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-sva:3.54.0--6dce402d4a7742a0':
        'community.wave.seqera.io/library/bioconductor-sva:3.54.0--6ac32dfb26a67399' }"

    input:
    val(condition)
    val(batch)
    tuple val(meta), path(derfinder_obj)

    output:
    tuple val(meta), path ("*.csv")  , emit: count_matrix
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'batchcorrection_script.R'

}
