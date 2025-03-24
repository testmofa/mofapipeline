process MOFA_PLOT {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bioconductor-mofa2=1.16 conda-forge::r-base=4.4.3 bioconda::bioconductor-mofadata=1.22 conda-forge::r-data.table=1.17.0 conda-forge::r-ggplot2=3.5.1"

    input:
    tuple val(meta), path(model)

    output:
    tuple val(meta), path("*.pdf"), emit: plots
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    library(MOFA2)
    library(ggplot2)

    # Load the MOFA model
    model <- load_model("${model}")

    # Plot variance explained
    pdf("${prefix}_variance_explained.pdf")
    plot_variance_explained(model, max_r2 = 15)
    dev.off()

    # Plot factor correlations
    pdf("${prefix}_factor_correlations.pdf")
    plot_factor_cor(model)
    dev.off()

    # Plot Top weights
    pdf("${prefix}_weights.pdf")
    plot_top_weights(model, view = "Drugs",
    factor 1, nfeatures = 10, scale = T)
    dev.off()

    # Plot factors
    pdf("${prefix}_factors.pdf")
    plot_factors(model, factors = 1:15)
    dev.off()

    # Create versions file
    writeLines(
        c(
            '"${task.process}":',
            paste0('    mofa: ', packageVersion("MOFA2")),
            paste0('    r-ggplot2: ', packageVersion("ggplot2"))
        ),
        "versions.yml"
    )
    """
}
