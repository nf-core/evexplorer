process DRE {
    // tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-edger_r-ggplot2_r-ggrepel_r-pheatmap_r-rcolorbrewer:6ab09c4ae74ac8f3':
        'community.wave.seqera.io/library/bioconductor-edger_r-ggplot2_r-ggrepel_r-pheatmap_r-rcolorbrewer:fe6f358940ad1878' }"

    input:
    val(condition)
    tuple val(meta), path(path_count_matrix)

    output:
    tuple val(meta), path ("*.csv"), emit: count_matrix
    tuple val(meta), path ("*.pdf"), emit: DE_analysis_plots
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'dre_script.R'


}
