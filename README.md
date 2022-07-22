[![Build Status](https://app.travis-ci.com/ncats/IntLIM.svg?branch=liz_dev)](https://app.travis-ci.com/ncats/IntLIM)

# IntLIM:  Integration through LInear Modeling

# IntLIM app is accessible via a server (no installation needed!).
Please [click here](https://intlim.ncats.io/).  And let us know if additional functionalities would be useful (see contact info below).

## IntLIM

Interpretation of metabolomics data is very challenging.  Yet it can be eased through integration of metabolomics with other ‘omics’ data. The IntLIM package, which includes a user-friendly RShiny web app, aims to integrate multiple types of omics data.  Unlike other approaches, IntLIM is focused on understanding how specific analyte associations are affected by phenotypic features.  To this end, we develop a linear modeling approach that describes how analyte associations are affected by phenotype.  The workflow involves the following steps: 1) input analyte level (e.g., expression or abundance) data files, 2) filter data sets by analyte level and imputed values, 3) run the linear model to extract FDR-adjusted interaction p-values, 4) filter results by p-values, interaction coefficient percentile, r-squared value, and Spearman correlation differences, and 5) plot/visualize specific analyte associations. 

An example data set is provided within the package, and is a subset of the NCI-60 gene expression and metabolomics data (https://wiki.nci.nih.gov/display/NCIDTPdata/Molecular+Target+Data).  The vignette outlines how to run the workflow. More details can be found in <a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2085-6" target="_blank"> our publication "IntLIM: integration using linear models of metabolomics and gene expression data"</a>.

## Citation
If you use IntLIM, please cite the following work:

Siddiqui JK, Baskin E, Liu M, Cantemir-Stone CZ, Zhang B, Bonneville R, McElroy JP, Coombes KR, Mathé EA. IntLIM: integration using linear models of metabolomics and gene expression data. BMC Bioinformatics. 2018 Mar 5;19(1):81. doi: 10.1186/s12859-018-2085-6.

PMID: 29506475; PMCID: PMC5838881 DOI: 10.1186/s12859-018-2085-6

To access, [click here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5838881/)


## IntLIM prerequisites

IntLIM is an R package and can be run on version >= 3.5.0. 

## Installation from Github

To install IntLIM, simply type the following in the R terminal:

```
install.packages("devtools")
library(devtools)
devtools::install_github("ncats/IntLIM")
```
## Vignette

A detailed vignette can be found [here](https://ncats.github.io/IntLIM/IntLIM2.0_Vignette.html).

## Formatted Data and Analysis Codes

Formatted data and codes to reproduce the NCI-60 and breast cancer analyses can be obtained from the following GitHub repository:

[https://github.com/ncats/IntLIMVignettes](https://github.com/ncats/IntLIMVignettes)


## Running IntLIM's user-friendly web app:

The package functions can be run directly in the R console.  
Alternatively, to launch the web app, type the following in your R console:

```
library(IntLIM)
runIntLIMApp()
```

## Contact

If you encounter any problems running on the software, or find installation problems or bugs, please start an issue on the Issues tab or email Ewy Mathe at Ewy.Mathe@nih.gov or Tara Eicher at Tara.Eicher@nih.gov.  We are also very open to any comments, including how we can improve and ameliorate the package.
