/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GUNZIP                      } from '../modules/nf-core/gunzip/main'
include { WINDOWMASKER_USTAT          } from '../modules/nf-core/windowmasker/ustat/main'
include { WINDOWMASKER_MKCOUNTS       } from '../modules/nf-core/windowmasker/mkcounts/main'
include { REPEATMODELER_REPEATMODELER } from '../modules/nf-core/repeatmodeler/repeatmodeler/main'
include { REPEATMODELER_MASKER as REPEATMODELER_MASKER_DFAM          } from '../modules/local/repeatmasker/main'
include { REPEATMODELER_MASKER as REPEATMODELER_MASKER_REPEATMODELER } from '../modules/local/repeatmasker/main'
include { REPEATMODELER_MASKER as REPEATMODELER_MASKER_EXTERNAL      } from '../modules/local/repeatmasker/main'
include { REPEATMODELER_BUILDDATABASE } from '../modules/nf-core/repeatmodeler/builddatabase/main'
include { TANTAN                      } from '../modules/local/tantan.nf'
include { BEDTOOLS_CUSTOM             } from '../modules/local/bedtools.nf'
include { CUSTOMMODULE                } from '../modules/local/custommodule.nf'
include { SEQTK_CUTN as TANTAN_BED        } from '../modules/local/seqtk.nf'
include { SEQTK_CUTN as WINDOWMASKER_BED  } from '../modules/local/seqtk.nf'
include { SEQTK_CUTN as REPEATMODELER_BED } from '../modules/local/seqtk.nf'
include { GFASTATS as GFSTTANTAN      } from '../modules/nf-core/gfastats/main'
include { GFASTATS as GFSTRMSK_DFAM   } from '../modules/nf-core/gfastats/main'
include { GFASTATS as GFSTRMSK_RMOD   } from '../modules/nf-core/gfastats/main'
include { GFASTATS as GFSTRMSK_EXTR   } from '../modules/nf-core/gfastats/main'
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

    if (params.gzipped_input) {
        GUNZIP(ch_samplesheet)
        input_genomes = GUNZIP.out.gunzip
    } else {
        input_genomes = ch_samplesheet
    }

    // Simple tandem repeat masking with tantan
    //
    TANTAN ( input_genomes )
    GFSTTANTAN ( TANTAN.out.masked_fa )
    TANTAN_BED ( TANTAN.out.masked_fa )

    // De novo repeat detection with WindowMasker
    //
    WINDOWMASKER_MKCOUNTS ( input_genomes )
    WINDOWMASKER_USTAT ( WINDOWMASKER_MKCOUNTS.out.counts.join(input_genomes) )
    GFSTWINDOWMASK ( WINDOWMASKER_USTAT.out.intervals )
    WINDOWMASKER_BED ( WINDOWMASKER_USTAT.out.intervals )

    // De novo repeat discovery and detection with RepeatModeller and RepeatMasker
    //
    REPEATMODELER_BUILDDATABASE ( input_genomes )
    REPEATMODELER_REPEATMODELER ( REPEATMODELER_BUILDDATABASE.out.db )
    REPEATMODELER_MASKER_REPEATMODELER (
        REPEATMODELER_REPEATMODELER.out.fasta
            .join(input_genomes)
            .map {meta, fasta, ref -> [ [id:"${meta.id}_REPM", id_old:meta.id] , fasta, ref ] },
        []
    )
    GFSTRMSK_RMOD ( REPEATMODELER_MASKER_REPEATMODELER.out.fasta )
    REPEATMODELER_BED ( REPEATMODELER_MASKER_REPEATMODELER.out.fasta.map {meta, fasta -> [ [id:meta.id_old], fasta ]} )

    // Repeat detection with DFAM and RepeatMasker
    //
    GFSTRMSK_DFAM_maybeout = channel.empty()
    if (params.taxon) {
        REPEATMODELER_MASKER_DFAM (
            input_genomes.map {meta, fasta -> [ [id:"${meta.id}_DFAM"] , [], fasta ] },
            params.taxon
        )
        GFSTRMSK_DFAM ( REPEATMODELER_MASKER_DFAM.out.fasta )
        GFSTRMSK_DFAM_maybeout = GFSTRMSK_DFAM.out.assembly_summary
    }

    GFSTRMSK_EXTR_maybeout = channel.empty()
    if (params.repeatlib) {
        REPEATMODELER_MASKER_EXTERNAL (
            input_genomes.map {meta, ref -> [ [id:"${meta.id}_EXTR"] , file(params.repeatlib, checkIfExists:true), ref ] },
            []
        )
        GFSTRMSK_EXTR ( REPEATMODELER_MASKER_EXTERNAL.out.fasta )
        GFSTRMSK_EXTR_maybeout = GFSTRMSK_EXTR.out.assembly_summary
    }

    //
    // MODULE: CUSTOMMODULE
    //
    CUSTOMMODULE ( channel.empty()
        .mix(    GFSTTANTAN.out.assembly_summary.map {it[1]} )
        .mix(GFSTWINDOWMASK.out.assembly_summary.map {it[1]} )
        .mix( GFSTRMSK_RMOD.out.assembly_summary.map {it[1]} )
        .mix(GFSTRMSK_DFAM_maybeout             .map {it[1]} )
        .mix(GFSTRMSK_EXTR_maybeout             .map {it[1]} )
        .collect()
    )
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOMMODULE.out.tsv)

    //
    // MODULE: bedtools_custom
    //
    BEDTOOLS_CUSTOM (
        input_genomes
            .join(TANTAN_BED.out.bed_gz)
            .join(WINDOWMASKER_BED.out.bed_gz)
            .join(REPEATMODELER_BED.out.bed_gz)
    )

    ch_versions = ch_versions
        .mix(WINDOWMASKER_MKCOUNTS.out.versions.first())
        .mix(TANTAN.out.versions.first())
        .mix(REPEATMODELER_REPEATMODELER.out.versions.first())
        .mix(GFSTWINDOWMASK.out.versions.first())
        .mix(TANTAN_BED.out.versions.first())
        .mix(BEDTOOLS_CUSTOM.out.versions.first())

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
