/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC as FASTQC_PRE_TRIM } from '../modules/nf-core/fastqc/main'
include { PYCHOPPER } from '../modules/nf-core/pychopper/main'
include { MINIMAP2_ALIGN } from '../modules/nf-core/minimap2/align/main'
include { MINIMAP2_INDEX } from '../modules/nf-core/minimap2/index/main'
include { FASTQC as FASTQC_POST_TRIM } from '../modules/nf-core/fastqc/main'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'
include { STAR_ALIGN                } from '../modules/nf-core/star/align/main'
include { SAMTOOLS_SORT             } from '../modules/nf-core/samtools/sort/main'      
include { SAMTOOLS_INDEX            } from '../modules/nf-core/samtools/index/main'     
include { DERFINDER                 } from '../modules/local/derfinder/main'
include { BATCHCORRECTION           } from '../modules/local/batchcorrection/main'      
include { DRE                       } from '../modules/local/dre/main'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_evexplorer_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow LONG_READS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta  // reference genome
    ch_gtf_1
    ch_gtf_2
    ch_chr_names 

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    


    //
    // MODULE: Index reference genome
    //
    MINIMAP2_INDEX ( ch_fasta )


    ch_minimap2_index = MINIMAP2_INDEX.out.index.collect()
    ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions.first())

    //
    // MODULE: Run FastQC before trimming
    //
    FASTQC_PRE_TRIM (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_PRE_TRIM.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_PRE_TRIM.out.versions.first())

    //
    // MODULE: Run pychopper for trimming
    //
    PYCHOPPER (
        ch_samplesheet
    )
 ch_trimmed_reads = PYCHOPPER.out.fastq
 ch_versions = ch_versions.mix(PYCHOPPER.out.versions.first())

    //
    // MODULE: Run FastQC after trimming
    //
    FASTQC_POST_TRIM (
       ch_trimmed_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST_TRIM.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_POST_TRIM.out.versions)


    //
    // MODULE: Aligning reads
    //
 MINIMAP2_ALIGN (
    ch_trimmed_reads,
    ch_minimap2_index,
     true,          // bam_format
     'bai',         // bam_index_extension
     false,         // cigar_paf_format
     true 
 )

// Collect all BAM outputs in a single channel
ch_sorted_bam = MINIMAP2_ALIGN.out.bam
ch_bam_index  = MINIMAP2_ALIGN.out.index


// grouping all bam files
ch_grouped_bam = ch_sorted_bam
    .map { meta, bam -> tuple([id: 'sample1', single_end: true], bam) }
    .groupTuple()

// grouping all indices together
ch_grouped_index = ch_bam_index
    .map { meta, index -> tuple([id: 'sample1', single_end: true], index) }
    .groupTuple()


    //
    // MODULE: generating counts from the bam files
    //
DERFINDER(
  ch_grouped_bam,
  ch_grouped_index,
  ch_chr_names,
  ch_gtf_1,
  ch_gtf_2
   )


ch_samplesheet
    .map { meta, _ -> tuple(meta.id, meta.sample_batch) }
    .collect()
    .set { ch_batch }

ch_derfinder = DERFINDER.out.Rda


    //
    // MODULE: applying batch correction 
    //
  BATCHCORRECTION (
      ch_batch,
      ch_derfinder
                 )
				 
ch_samplesheet
    .map { meta, _ -> tuple(meta.id, meta.sample_batch, meta.sample_cond ) }
    .collect()
    .set { ch_condition }

ch_count_matrix= BATCHCORRECTION.out.count_matrix


    //
    // MODULE: performing differential RNA expression 
    //
DRE (
ch_condition,
ch_count_matrix
  )
  
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'evexplorer_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

