library(devtools)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MultiDataSet", update = TRUE)

install_github("ncats/IntLIM", force = TRUE)
library(IntLIM)
