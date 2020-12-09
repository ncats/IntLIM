library(devtools)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MultiDataSet", update = TRUE)

install.packages("shinyFiles")
