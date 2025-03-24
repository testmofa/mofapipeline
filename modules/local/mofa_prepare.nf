process MOFA_PREPARE {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bioconductor-mofa2=1.16 conda-forge::r-base=4.4.3 bioconda::bioconductor-mofadata=1.22 conda-forge::r-data.table=1.17.0 conda-forge::r-ggplot2=3.5.1"

    input:
    tuple val(meta), path(omics_files)

    output:
    tuple val(meta), path("*.rds"), emit: prepared_data
    path "versions.yml"           , emit: versions

    script:
    """
    #!/usr/bin/env Rscript
    library(MOFA2)
    library(MOFAdata)
    library(data.table)
    library(ggplot2)
    library(tidyverse)

    # Your data preparation code here
    
    utils::data("CLL_data")
    lapply(CLL_data,dim)
    
    CLL_metadata <- fread("ftp://ftp.ebi.ac.uk/pub/databases/mofa/cll_vignette/sample_metadata.txt")
    
    # Save prepared data
    saveRDS(CLL_data, CLL_metadata file = "${meta.id}_prepared_data.rds")

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
