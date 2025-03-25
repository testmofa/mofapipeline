#!/usr/bin/env nextflow

//Channel.fromPath('C:\Users\faria\OneDrive\Documentos\mofapipeline\test_mofa\*.csv') // Matches all .txt files in the "data" folder
  //  .set { file_list }

process MOFA_PREPARE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bioconductor-mofa2=1.16 conda-forge::r-base=4.4.3 bioconda::bioconductor-mofadata=1.22 conda-forge::r-data.table=1.17.0 conda-forge::r-ggplot2=3.5.1"

    input:
    //file input_files from params.file_list 
    tuple val(meta_drugs), path(drugs_data)
    tuple val(meta_meth), path(meth_data)
    tuple val(meta_mRNA), path(mRNA_data)
    tuple val(meta_muta), path(muta_data)


    output:
    tuple val(meta), path("*.rds"), emit: prepared_data
    path "versions.yml"           , emit: versions

    script:
    """
    #!/usr/bin/env Rscript ${input_files}

    library(MOFA2)
    library(MOFAdata)
    library(data.table)
    library(ggplot2)
    library(tidyverse)
    
    table_drug <- read.table("$meta_drug")
    table_meth <- read.table("$meta_meth")
    table_mRNA <- read.table("$meta_mRNA")
    table_muta <- read.table("$meta_muta")

    MOFA_LIST <- list()    
    MOFA_LIST["mRNA"] <- table_drug
    MOFA_LIST["Mutations"] <- table_meth    
    MOFA_LIST["Methylation"] <- table_mRNA
    MOFA_LIST["Drugs"] <- table_muta

    # Save prepared data
    saveRDS(MOFA_LIST file = "prepared_data.rds")

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

workflow { 
    prepare_ch = MOFA_PREPARE(
        params.drugs_data,
        params.meth_data,
        params.mRNA_data,
        params.muta_data
    ) 
    run_ch = MOFA_RUN(MOFA_PREPARE.out.prepared_data)
    plot_ch = MOFA_PLOT(MOFA_RUN.out.model) 
    plot_ch.view { it } 
}