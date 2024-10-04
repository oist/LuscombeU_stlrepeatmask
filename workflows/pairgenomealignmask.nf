/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GUNZIP                                             } from '../modules/nf-core/gunzip/main'

include { TANTAN               as TANTAN_MASK                } from '../modules/local/tantan.nf'
include { GFASTATS             as TANTAN_STATS               } from '../modules/nf-core/gfastats/main'
include { SEQTK_CUTN           as TANTAN_BED                 } from '../modules/local/seqtk.nf'

include { WINDOWMASKER_MASK                                  } from '../modules/local/windowmasker.nf'
include { GFASTATS             as WINDOWMASKER_STATS         } from '../modules/nf-core/gfastats/main'
include { SEQTK_CUTN           as WINDOWMASKER_BED           } from '../modules/local/seqtk.nf'

include { REPEATMODELER_REPEATMODELER                        } from '../modules/nf-core/repeatmodeler/repeatmodeler/main'
include { REPEATMODELER_BUILDDATABASE                        } from '../modules/nf-core/repeatmodeler/builddatabase/main'
include { REPEATMODELER_MASKER as REPEATMODELER_REPEATMASKER } from '../modules/local/repeatmasker.nf'
include { GFASTATS             as REPEATMODELER_STATS        } from '../modules/nf-core/gfastats/main'
include { SEQTK_CUTN           as REPEATMODELER_BED          } from '../modules/local/seqtk.nf'

include { REPEATMODELER_MASKER as DFAM_REPEATMASKER          } from '../modules/local/repeatmasker.nf'
include { GFASTATS             as DFAM_STATS                 } from '../modules/nf-core/gfastats/main'
include { SEQTK_CUTN           as DFAM_BED                   } from '../modules/local/seqtk.nf'

include { REPEATMODELER_MASKER as EXTLIB_REPEATMASKER        } from '../modules/local/repeatmasker.nf'
include { GFASTATS             as EXTLIB_STATS               } from '../modules/nf-core/gfastats/main'
include { SEQTK_CUTN           as EXTLIB_BED                 } from '../modules/local/seqtk.nf'

include { MERGE_REPM_RESULTS   as MERGEDMASKS_REPM           } from '../modules/local/merge_repeatmasker_results.nf'
include { GFASTATS             as MERGEDMASKS_REPM_STATS     } from '../modules/nf-core/gfastats/main'

include { MERGE_MASKS          as MERGEDMASKS_ALL            } from '../modules/local/mergemasks.nf'
include { GFASTATS             as MERGEDMASKS_ALL_STATS      } from '../modules/nf-core/gfastats/main'

include { MULTIQC_SOFTMASK_STATS                             } from '../modules/local/multiqc_softmask_statistics.nf'
include { MULTIQC_SOFTMASK_OVERLAPS                          } from '../modules/local/multiqc_softmask_overlaps.nf'

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
    TANTAN_MASK  ( input_genomes.map { meta, ref -> [ [id:"${meta.id}_tantan", key:meta.id], ref ] } )
    TANTAN_STATS ( TANTAN_MASK.out.masked_fa )
    TANTAN_BED   ( TANTAN_MASK.out.masked_fa )

    // De novo repeat detection with WindowMasker
    //
    WINDOWMASKER_MASK  ( input_genomes.map { meta, ref -> [ [id:"${meta.id}_windowmasker", key:meta.id], ref ] } )
    WINDOWMASKER_STATS ( WINDOWMASKER_MASK.out.masked_fa )
    WINDOWMASKER_BED   ( WINDOWMASKER_MASK.out.masked_fa )

    // De novo repeat discovery and detection with RepeatModeller and RepeatMasker
    //
    REPEATMODELER_BUILDDATABASE ( input_genomes )
    REPEATMODELER_REPEATMODELER ( REPEATMODELER_BUILDDATABASE.out.db )
    REPEATMODELER_REPEATMASKER  (
        REPEATMODELER_REPEATMODELER.out.fasta
            .join(input_genomes)
            .map {meta, fasta, ref -> [ [id:"${meta.id}_REPM", key:meta.id] , fasta, ref ] },
        []
    )
    REPEATMODELER_STATS         ( REPEATMODELER_REPEATMASKER.out.fasta )
    REPEATMODELER_BED           ( REPEATMODELER_REPEATMASKER.out.fasta )

    // Repeat detection with DFAM and RepeatMasker
    //
    DFAM_BED_maybeout  = channel.empty()
    DFAM_STATS_maybeout = channel.empty()
    if (params.taxon) {
        DFAM_REPEATMASKER (
            input_genomes.map {meta, fasta -> [ [id:"${meta.id}_DFAM", key:meta.id] , [], fasta ] },
            params.taxon
        )
        DFAM_STATS ( DFAM_REPEATMASKER.out.fasta )
        DFAM_BED   ( DFAM_REPEATMASKER.out.fasta )
        DFAM_BED_maybeout   = DFAM_BED.out
        DFAM_STATS_maybeout = DFAM_STATS.out.assembly_summary
    }

    // Repeat detection with RepeatMasker and an external library of repeats
    //
    EXTLIB_REPM_maybeout  = channel.empty()
    EXTLIB_STATS_maybeout = channel.empty()
    if (params.repeatlib) {
        EXTLIB_REPEATMASKER (
            input_genomes.map {meta, ref -> [ [id:"${meta.id}_EXTR", key:meta.id] , file(params.repeatlib, checkIfExists:true), ref ] },
            []
        )
        EXTLIB_STATS ( EXTLIB_REPEATMASKER.out.fasta )
        EXTLIB_BED   ( EXTLIB_REPEATMASKER.out.fasta )
        EXTLIB_BED_maybeout   = EXTLIB_BED.out
        EXTLIB_STATS_maybeout = EXTLIB_STATS.out.assembly_summary
    }

    // Comparing and merging soft masks from each repeatmasker run
    //
    MERGEDMASKS_REPM (
        input_genomes
            .join(REPEATMODELER_BED.out.bed_gz.map {meta, bed -> [ [id:meta.key ] , bed ] } )
            .join(    DFAM_BED_maybeout.bed_gz.map {meta, bed -> [ [id:meta.key ] , bed ] } )
            .join(  EXTLIB_BED_maybeout.bed_gz.map {meta, bed -> [ [id:meta.key ] , bed ] } )
    )
    MERGEDMASKS_REPM_STATS ( MERGEDMASKS_REPM.out.fasta )

    // Comparing and merging soft masks from each software.
    //
    MERGEDMASKS_ALL (
        input_genomes
            .join(TANTAN_BED.out.bed_gz.map       {meta, bed -> [ [id:meta.key ] , bed ] } )
            .join(WINDOWMASKER_BED.out.bed_gz.map {meta, bed -> [ [id:meta.key ] , bed ] } )
            .join(MERGEDMASKS_REPM.out.repm_all_bed_gz                                     )
    )
    MERGEDMASKS_ALL_STATS ( MERGEDMASKS_ALL.out.fasta )

    // Aggregation of statistics
    //
    MULTIQC_SOFTMASK_STATS ( channel.empty()
        .mix(          TANTAN_STATS.out.assembly_summary.map {it[1]} )
        .mix(    WINDOWMASKER_STATS.out.assembly_summary.map {it[1]} )
        .mix(   REPEATMODELER_STATS.out.assembly_summary.map {it[1]} )
        .mix(            DFAM_STATS_maybeout            .map {it[1]} )
        .mix(          EXTLIB_STATS_maybeout            .map {it[1]} )
        .mix( MERGEDMASKS_REPM_STATS.out.assembly_summary.map{it[1]} )
        .mix(  MERGEDMASKS_ALL_STATS.out.assembly_summary.map{it[1]} )
        .collect()
    )
    ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_SOFTMASK_STATS.out.tsv)

    // Aggregation of statistics (Jaccard indices)
    //
    MULTIQC_SOFTMASK_OVERLAPS ( MERGEDMASKS_ALL.out.txt.map{it[1]}.collect() )
    ch_multiqc_files = ch_multiqc_files.mix(MULTIQC_SOFTMASK_OVERLAPS.out.tsv)


    // Collect software versions
    //
    ch_versions = ch_versions
        .mix(WINDOWMASKER_MASK.out.versions.first())
        .mix(TANTAN_MASK.out.versions.first())
        .mix(REPEATMODELER_REPEATMODELER.out.versions.first())
        .mix(WINDOWMASKER_STATS.out.versions.first())
        .mix(TANTAN_BED.out.versions.first())
        .mix(MERGEDMASKS_ALL.out.versions.first())

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    // MultiQC
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
