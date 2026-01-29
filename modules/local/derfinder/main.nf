process DERFINDER {
    // tag '$bam'
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-derfinder_r-caret_r-data.table_r-devtools_pruned:db22279c2b4717b0':
        'community.wave.seqera.io/library/bioconductor-derfinder_r-caret_r-data.table_r-devtools_pruned:570fa1e65dc6e12d' }"

    input:
    tuple val(meta), path(path_bam)
    tuple val(meta), path(path_bai)
    path chr_names
    path gtf_1
    path gtf_2

    output:
    tuple val(meta), path ("*.Rda")  , emit: Rda 
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    template 'derfinder_script.R'    
    

}
