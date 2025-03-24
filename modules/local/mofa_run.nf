process MOFA_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bioconductor-mofa2=1.16 conda-forge::r-base=4.4.3 bioconda::bioconductor-mofadata=1.22 conda-forge::r-data.table=1.17.0 conda-forge::r-ggplot2=3.5.1"

    input:
    tuple val(meta), path(prepared_data)

    output:
    tuple val(meta), path("*.rds") , emit: model
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    library(MOFA2)

    # Load prepared data
    data <- readRDS("${prepared_data}")

    # Create and train the MOFA model
    model <- create_mofa(data)
    model <- prepare_mofa(model)
    model <- run_mofa(model)

    # Save the model
    saveRDS(model, file = "${prefix}_mofa_model.rds")

    # Create versions file
    writeLines(
        c(
            '"${task.process}":',
            paste0('    mofa: ', packageVersion("MOFA2"))
        ),
        "versions.yml"
    )
    """
}
