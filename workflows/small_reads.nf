/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC as FASTQC_PRE_TRIM } from '../modules/nf-core/fastqc/main'
include { CUTADAPT                  } from '../modules/nf-core/cutadapt/main'
include { FASTQC as FASTQC_POST_TRIM } from '../modules/nf-core/fastqc/main'
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'
include { STAR_GENOMEGENERATE       } from '../modules/nf-core/star/genomegenerate/main'
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

workflow SMALL_READS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta
    ch_gtf_1
    ch_gtf_2
    ch_chr_names

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
       
       //
       // MODULE: STAR_GENOMEGENERATE
       //
  STAR_GENOMEGENERATE (
       ch_fasta,
         [ [:], [] ]
)


ch_star_index = STAR_GENOMEGENERATE.out.index.collect()
ch_versions = ch_versions.mix(STAR_GENOMEGENERATE.out.versions.first())

    //
    // MODULE: Run FastQC before trimming
    //
    FASTQC_PRE_TRIM (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_PRE_TRIM.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_PRE_TRIM.out.versions.first())

    //
    // MODULE: Run Cutadapt for trimming
    //
    CUTADAPT (
        ch_samplesheet
    )
    ch_trimmed_reads = CUTADAPT.out.reads
    ch_versions = ch_versions.mix(CUTADAPT.out.versions.first())

    //
    // MODULE: Run FastQC after trimming
    //
    FASTQC_POST_TRIM (
       ch_trimmed_reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_POST_TRIM.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC_POST_TRIM.out.versions)


    //
    // MODULE: STAR_ALIGN
    //
    STAR_ALIGN (
       ch_trimmed_reads,
       ch_star_index,
       [ [:], [] ],
       [ [:], [] ],
       [ [:], [] ],
       [ [:], [] ]
    )

    // Collect all BAM outputs in a single channel
    ch_aligned_bam = STAR_ALIGN.out.bam

   
    // MODULE: Samtools Sort
    SAMTOOLS_SORT (
        ch_aligned_bam,
        ch_fasta
    )
    ch_sorted_bam = SAMTOOLS_SORT.out.bam
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())



    // MODULE: Samtools Index
    SAMTOOLS_INDEX (
        ch_sorted_bam
    )
    ch_bam_index = SAMTOOLS_INDEX.out.bai
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())



     ch_grouped_bam = ch_sorted_bam
    .map { meta, bam -> tuple([id: 'sample1', single_end: true], bam) }
    .groupTuple()

     ch_grouped_index = ch_bam_index
    .map { meta, index -> tuple([id: 'sample1', single_end: true], index) }
    .groupTuple()

     //
     // MODULE: DERFINDER to generate counts from bam files
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


    ch_samplesheet
    .map { meta, _ -> tuple(meta.id, meta.sample_batch, meta.sample_cond ) }
    .collect()
    .set { ch_condition }

    ch_derfinder = DERFINDER.out.Rda

     
     //
     // MODULE: BATCHCORRECTION
     //
     BATCHCORRECTION (
      ch_condition,
      ch_batch,
      ch_derfinder
                 )


   ch_count_matrix= BATCHCORRECTION.out.count_matrix

    //
    // MODULE: performing differential RNA expression
    //
    DRE (
    ch_condition,
    ch_count_matrix
    ) 
   
 
 
 
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

