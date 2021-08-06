#' MetaboliteSet class extenstion for MultiDataSet object
#'
#' This class is specific to metabolomics data and will then be stored into the class
#' \code{MultiDataSet}, which stores multiple dataset types (e.g. metabolite and gene levels
#' in this case.
#'
#' @name MetaboliteSet-class
#' @rdname MetaboliteSet-class
#' @exportClass MetaboliteSet
#' @slot eSet List of eSet elements

if(!require("Biobase")) install.packages("Biobase")

methods::setClass (
	Class = "MetaboliteSet",
	contains = "eSet",
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
#' @slot warnings a message of whether genes and metabolites have 0 standard deviation
#' @slot stype column name that represents sample type (by default, it will be used
#' in the interaction term). Only 2 categories are currently supported.
#' @slot outcome outcome is either 'metabolite' or 'gene'
#' @slot independent.var.type independent variable type (either 'metabolite or 'gene')
#' @slot covar describing additional variables and the class they form
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
		warnings="character",
		stype="character",
		outcome="character",
		independent.var.type="character",
		covar="data.frame")
	)
#' ModelInput class
#'
#' @name ModelInput-class
#' @rdname ModelInput-class
#' @exportClass ModelInput
#' @slot A.hat The Laplacian of a line graph built from the co-regulation graphs, where 
#' each node corresponds to a pair of analytes.
#' @slot node.wise.prediction graph of co-regulation results from IntLIM
#' @slot true.phenotypes data frame of true phenotypes given analyte pairs
#' and clinical covariates
#' @slot outcome.type "numeric" or "categorical"
methods::setClass(
  Class="ModelInput",
  representation(A.hat="matrix",
                 node.wise.prediction="matrix",
                 true.phenotypes="numeric",
                 outcome.type="character")
)
#' PoolingFilter class
#'
#' @name PoolingFilter-class
#' @rdname PoolingFilter-class
#' @exportClass PoolingFilter
#' @slot filter The filter that maps each dimension of the input to one of k clusters.
#' @slot filter.type One of "mean", "min", "median", or "max".
#' @slot cluster.sizes A vector of cluster sizes (by number of dimensions mapping
#' @slot individual.filters A list of filters corresponding to each individual sample.
#' to the cluster)
methods::setClass(
  Class="PoolingFilter",
  representation(filter="matrix",
                 filter.type="character",
                 cluster.sizes="matrix",
                 individual.filters="list")
)
#' ModelResults class
#'
#' @name ModelResults-class
#' @rdname ModelResults-class
#' @exportClass ModelResults
#' @slot model.input An object of the modelInput class.
#' @slot pooling.filter A matrix that pools convolution results.
#' @slot iteration.tracking A data frame to track the iteration, weight, and error
#' values for each iteration of training.
#' @slot max.iterations Maximum number of iterations.
#' @slot convergence.cutoff Cutoff for convergence.
#' @slot learning.rate Learning rate used during training.
#' @slot activation.type Character value. Must be "sigmoid", "tanh", or "softmax".
#' @slot current.weights Weights used in the current iteration.
#' @slot previous.weights Weights used in the previous iteration.
#' @slot current.gradient Gradient calculated for this iteration.
#' @slot weights.after.pooling Whether to include the weights after the pooling
#' layer (as opposed to before). Must be a boolean.
#' @slot outcome.prediction The prediction of the outcome.
#' @slot optimization.type One of "BGD", "momentum", "adagrad", or "adam"
#' @slot current.iteration Iteration (changes at each time step)
#' @slot previous.momentum Momentum value used in momentum optimization
#' @slot previous.update.vector Previous value used to update weights, used in ADAM
#' optimization
#' @slot sum.square.gradients Sum of squared gradients over iterations, used in
#' Adagrad optimization
methods::setClass(
  Class="ModelResults",
  representation(model.input="ModelInput",
                 pooling.filter="PoolingFilter",
                 iteration.tracking="data.frame",
                 max.iterations="numeric",
                 convergence.cutoff="numeric",
                 learning.rate="numeric",
                 activation.type="character",
                 previous.weights="matrix",
                 current.weights="matrix",
                 current.gradient="matrix",
                 previous.update.vector="matrix",
                 current.iteration="numeric",
                 weights.after.pooling="logical",
                 outcome.prediction="numeric",
                 optimization.type="character",
                 sum.square.gradients="matrix",
                 previous.momentum="matrix")
)


