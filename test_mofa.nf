#!/usr/bin/env nextflow

//Channel.fromPath('C:\Users\faria\OneDrive\Documentos\mofapipeline\test_mofa\*.csv') // Matches all .txt files in the "data" folder
  //  .set { file_list }

process MOFA_PREPARE {
    label 'process_medium'

    conda "conda-forge::r-base=4.4.1 bioconda::bioconductor-mofa2=1.16 bioconda::bioconductor-mofadata=1.22 conda-forge::r-data.table=1.17.0 conda-forge::r-ggplot2=3.5.1 conda-forge::r-tidyverse conda-forge::icu"

    input:
    //file input_files from params.file_list
    path(drugs_data)
    path(meth_data)
    path(mRNA_data)
    path(muta_data)


    output:
    path("*.rds"), emit: prepared_data
    path "versions.yml"           , emit: versions

    script:
    """
    #!/usr/bin/env Rscript

    library("MOFA2")
    library("MOFAdata")
    library("data.table")
    library("ggplot2")
    library("tidyverse")

    table_drug <- read.table("$drugs_data")
    table_meth <- read.table("$meth_data")
    table_mRNA <- read.table("$mRNA_data")
    table_muta <- read.table("$muta_data")

    MOFA_LIST <- list()
    MOFA_LIST["mRNA"] <- table_drug
    MOFA_LIST["Mutations"] <- table_meth
    MOFA_LIST["Methylation"] <- table_mRNA
    MOFA_LIST["Drugs"] <- table_muta

    # Save prepared data
    saveRDS(MOFA_LIST, file = "prepared_data.rds")

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

process MOFA_RUN {
    label 'process_medium'

    conda "bioconda::bioconductor-mofa2=1.16 conda-forge::r-base=4.4.3 bioconda::bioconductor-mofadata=1.22 conda-forge::r-data.table=1.17.0 conda-forge::r-ggplot2=3.5.1"

    input:
    path(prepared_data)

    output:
    path("*.rds") , emit: model
    path "versions.yml"            , emit: versions

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    library(MOFA2)

    # Load prepared data
    data <- readRDS("$prepared_data")

    # Create and train the MOFA model
    model <- create_mofa(data)
    model <- prepare_mofa(model)
    model <- run_mofa(model)

    # Save the model
    saveRDS(model, file = "mofa_model.rds")

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

process MOFA_PLOT {
    label 'process_medium'

    conda "bioconda::bioconductor-mofa2=1.16 conda-forge::r-base=4.4.3 bioconda::bioconductor-mofadata=1.22 conda-forge::r-data.table=1.17.0 conda-forge::r-ggplot2=3.5.1"

    input:
    path(model)

    output:
    path("*.pdf"), emit: plots
    path "versions.yml"           , emit: versions

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    """
    #!/usr/bin/env Rscript
    library(MOFA2)
    library(ggplot2)

    # Load the MOFA model
    model <- load_model("${model}")

    # Plot variance explained
    pdf("variance_explained.pdf")
    plot_variance_explained(model, max_r2 = 15)
    dev.off()

    # Plot factor correlations
    pdf("factor_correlations.pdf")
    plot_factor_cor(model)
    dev.off()

    # Plot Top weights
    pdf("weights.pdf")
    plot_top_weights(model, view = "Drugs",
    factor 1, nfeatures = 10, scale = T)
    dev.off()

    # Plot factors
    pdf("factors.pdf")
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

workflow {
    prepare_ch = MOFA_PREPARE(
        file(params.drugs_data),
        file(params.meth_data),
        file(params.mRNA_data),
        file(params.muta_data)
    )
    run_ch = MOFA_RUN(MOFA_PREPARE.out.prepared_data)
    plot_ch = MOFA_PLOT(MOFA_RUN.out.model)
}
