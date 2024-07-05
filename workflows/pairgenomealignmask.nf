/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { WINDOWMASKER_USTAT          } from '../modules/nf-core/windowmasker/ustat/main'
include { WINDOWMASKER_MKCOUNTS       } from '../modules/nf-core/windowmasker/mkcounts/main'
include { REPEATMODELER_REPEATMODELER } from '../modules/nf-core/repeatmodeler/repeatmodeler/main'
include { REPEATMODELER_MASKER        } from '../modules/nf-core/repeatmodeler/repeatmasker/main'
include { REPEATMODELER_BUILDDATABASE } from '../modules/nf-core/repeatmodeler/builddatabase/main'
include { TANTAN                      } from '../modules/local/tantan.nf'
include { GFASTATS as GFSTTANTAN      } from '../modules/nf-core/gfastats/main'
include { GFASTATS as GFSTREPEATMOD   } from '../modules/nf-core/gfastats/main'
include { GFASTATS as GFSTWINDOWMASK  } from '../modules/nf-core/gfastats/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap            } from 'plugin/nf-validation'
include { paramsSummaryMultiqc        } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText      } from '../subworkflows/local/utils_nfcore_pairgenomealignmask_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PAIRGENOMEALIGNMASK {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    //
    // MODULE: tantan
    //
    TANTAN (
        ch_samplesheet
    )
    
    //
    // MODULE: gfastats_tantan
    //
    GFSTTANTAN (
        TANTAN.out.masked_fa
    )
    
    // MODULE: repeatmodeler_builddatabase
    //
    REPEATMODELER_BUILDDATABASE (
        ch_samplesheet
    )

    //
    // MODULE: repeatmodeler_repeatmodeler
    //
    REPEATMODELER_REPEATMODELER (
        REPEATMODELER_BUILDDATABASE.out.db
    )

    //
    // MODULE: repeatmodeler_repeatmasker
    //
    REPEATMODELER_MASKER (
        REPEATMODELER_REPEATMODELER.out.fasta,
        ch_samplesheet
    )

    //
    // MODULE: gfastats_repeatmodeler
    //
    GFSTREPEATMOD (
        REPEATMODELER_REPEATMODELER.out.fasta
    )

    //
    // MODULE: windowmasker_mkcounts
    //
    WINDOWMASKER_MKCOUNTS (
        ch_samplesheet
    )

    //
    // MODULE: windowmasker_ustat
    //
    WINDOWMASKER_USTAT (
        WINDOWMASKER_MKCOUNTS.out.counts.join(ch_samplesheet)
    )

    //
    // MODULE: gfastats_windowmasker
    //
    GFSTWINDOWMASK (
        WINDOWMASKER_USTAT.out.intervals
    )

    ch_multiqc_files = ch_multiqc_files.mix(WINDOWMASKER_MKCOUNTS.out.counts.collect{it[1]})
    ch_versions = ch_versions.mix(WINDOWMASKER_MKCOUNTS.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
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

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
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
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
