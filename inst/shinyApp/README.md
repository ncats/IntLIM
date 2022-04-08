# Welcome to the IntLIM shiny app!

The goal of this app is to provide users with a user-friendly platform for integrating multi-omics data.  Specifically, the software finds analyte relationships that are specific to a given phenotype (e.g. cancer vs non-cancer). For example, a given analyte pair could show a strong correlation in one phenotype (e.g. cancer) and no correlation in the other (e.g. non-cancer). 

More details can be found in <a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2085-6" target="_blank"> our publication "IntLIM: integration using linear models of metabolomics and gene expression data"</a>.

## Getting started (loading in data)

__*Please be sure that all files noted in the CSV file, including the CSV file, are in the same folder. Do not include path names in the filenames.*__

Users will need to input files for analyte levels for analyte type 1 (e.g. metabolite abundance data), analyte levels for analyte type 2 (e.g. gene expression data), sample meta-data, analyte type 1 meta-data (optional) and analyte type 2 meta-data (optional).  

Users also need to input a CSV file named 'input.csv' with two required columns: 'type' and 'filenames'.

The CSV file is expected to have the following 2 columns and 6 rows:

1. type,filenames
2. analyteType1,myfilename
3. analyteType2,myfilename
4. analyteType1MetaData,myfilename (optional)
5. analyteType2MetaData,myfilename (optional)
6. sampleMetaData,myfilename"

*NOTE*: For the ShinyApp, the meta-file must be named 'input.csv'

Note also that the input data files should be in a specific format:
- analyteType1: rows are analytes, columns are samples; the first row is assumed to have sample ids and these ids should be unique; the first column is assumed to have feature ids and those should be unique.
- analyteType2: rows are analytes, columns are samples; the first row is assumed to have sample ids and these ids should be unique; the first column is assumed to have feature ids and those should be unique.
- analyteType1MetaData: rows are analytes, features are columns
- analyteType2MetaData: rows are analytes, features are columns
- sampleMetaData: rows are samples, features are columns

*NOTE*: The first column of the sampleMetaData file is assumed to be the sample id, and those sample ids should match the *first row* of analyteType1 and analyteType2 (e.g. it is required that all sample ids in the analyteType1 and analyteType2 are also in the sampleMetaDatafile).

## Test data
The package includes a reduced set of the original NCI-60 dataset.  The CSV input file location for this test dataset can be located by typing the following in the R console:
```
     dir <- system.file("extdata", package="IntLIM", mustWork=TRUE)
     csvfile <- file.path(dir, "NCItestinput.csv")
     csvfile
```
Please see the vignette at [https://github.com/ncats/IntLIM/tree/liz_dev/vignettes/IntLimVignette.Rmd) for additional information.

In addition, additional NCI-60 and breast cancer demo datasets can be found at [https://github.com/ncats/IntLIM2.0ExtraDataVignettes](https://github.com/ncats/IntLIM2.0ExtraDataVignettes).

## Contact

If you have any questions, comments, or concerns on how to use IntLIM please contact Ewy Mathe at ewy.mathe@nih.gov or  Tara Eicher at tara.eicher@nih.gov.
