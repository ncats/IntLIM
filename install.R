library(devtools)

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("MultiDataSet", update = TRUE)

update.packages(ask = FALSE)

install_github("ncats/IntLIM", force = TRUE)
library(IntLIM)
