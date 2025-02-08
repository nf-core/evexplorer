process DGE {
    // tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-edger:4.4.0--dc3b5f57c26c3a8e':
        'community.wave.seqera.io/library/bioconductor-edger:4.4.0--fb6573efa683e367' }"

    input:
    tuple val(meta), path(path_sample_sheet)
    tuple val(meta), path(path_count_matrix)

    output:
    tuple val(meta), path ("*.csv"), emit: count_matrix
    tuple val(meta), path ("*.pdf"), emit: DE_analysis_plots
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'differentialgeneexpression.R'


}
