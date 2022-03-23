#' IntLimData class
#'
#' @name IntLimData-class
#' @rdname IntLimData-class
#' @exportClass IntLimData
#' @slot analyteType1 A matrix of abundance, expression, or other levels for
#' a specific type of analyte (e.g. protein abundance, metabolite abundance, or
#' gene expression)
#' @slot analyteType2 A second matrix of abundance, expression, or other levels for
#' a specific type of analyte (e.g. protein abundance, metabolite abundance, or
#' gene expression)
#' @slot analyteType1MetaData A data frame of metadata for analyte type 1.
#' @slot analyteType2MetaData A data frame of metadata for analyte type 2.
#' @slot sampleMetaData A data frame of covariate values from the patient data.
methods::setClass(
  Class="IntLimData",
  representation(analyteType1="matrix",
                 analyteType2="matrix",
                 analyteType1MetaData = "data.frame",
                 analyteType2MetaData = "data.frame",
                 sampleMetaData = "data.frame")
                 
)

#' IntLimResults class
#'
#' @name IntLimResults-class
#' @rdname IntLimResults-class
#' @exportClass IntLimResults
#' @slot interaction.pvalues matrix of interaction p-values
#' @slot interaction.adj.pvalues matrix of adjusted interaction pvalues
#' @slot interaction.coefficients matrix of interaction coefficients
#' @slot covariate.coefficients data frame of coefficients for each covariate
#' @slot covariate.pvalues data frame of p-values for each covariate
#' @slot model.rsquared matrix of r-squared values
#' @slot corr matrix of correlations in group 1 and 2
#' @slot filt.results data frame of filtered results
#' @slot warnings a message of whether analytes have 0 standard deviation
#' @slot stype column name that represents sample type (by default, it will be used
#' in the interaction term). Only 2 categories are currently supported.
#' @slot outcome outcome is either '1' or '2'
#' @slot independent.var.type independent variable type (either '1 or '2')
#' @slot covar describing additional variables and the class they form
#' @slot continuous "1" if outcome is continuous, "0" if not
methods::setClass(
	Class="IntLimResults",
	representation(interaction.pvalues="matrix",
		interaction.adj.pvalues="matrix",
		interaction.coefficients="matrix",
		model.rsquared = "matrix",
		covariate.pvalues = "data.frame",
		covariate.coefficients = "data.frame",
		corr="data.frame",
		filt.results="data.frame",
		warnings="list",
		stype="character",
		outcome="numeric",
		independent.var.type="numeric",
		covar="character",
		continuous="numeric")
	)