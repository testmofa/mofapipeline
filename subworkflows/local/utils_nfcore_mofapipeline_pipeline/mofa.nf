//
// MOFA (Multi-Omics Factor Analysis) subworkflow
//

include { MOFA_PREPARE } from '../../modules/local/mofa_prepare'
include { MOFA_RUN     } from '../../modules/local/mofa_run'
include { MOFA_PLOT    } from '../../modules/local/mofa_plot'

workflow MOFA {
    take:
    omics_data // channel: [ val(meta), [ omics_files ] ]

    main:
    ch_versions = Channel.empty()

    //
    // MODULE: Prepare data for MOFA
    //
    MOFA_PREPARE ( omics_data )
    ch_versions = ch_versions.mix(MOFA_PREPARE.out.versions)

    //
    // MODULE: Run MOFA
    //
    MOFA_RUN ( MOFA_PREPARE.out.prepared_data )
    ch_versions = ch_versions.mix(MOFA_RUN.out.versions)

    //
    // MODULE: Plot MOFA results
    //
    MOFA_PLOT ( MOFA_RUN.out.model )
    ch_versions = ch_versions.mix(MOFA_PLOT.out.versions)

    emit:
    model    = MOFA_RUN.out.model    // channel: [ val(meta), path(model) ]
    plots    = MOFA_PLOT.out.plots   // channel: [ val(meta), path(plots) ]
    versions = ch_versions           // channel: [ versions.yml ]
}
