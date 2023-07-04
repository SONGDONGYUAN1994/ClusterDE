#' Construct the synthetic null data
#'
#' \code{constructNull} takes the target data as the input and returns the corresponding synthetic null data.
#'
#' This function constructs the synthetic null data based on the target data (real data). The input is a expression matrix (gene by cell); the user should specify a distribution, which is usually Negative Binomial for count matrix.
#'
#' @param mat An expression matrix (gene by cell). It can be a regular matrix or a \code{sparseMatrix}.
#' @param family A string or a vector of strings of the distribution of your data.
#' Must be one of 'nb', 'binomial', 'poisson', 'zip', 'zinb' or 'gaussian', which represent 'poisson distribution',
#' 'negative binomial distribution', 'zero-inflated poisson distribution', 'zero-inflated negative binomail distribution',
#' and 'gaussian distribution' respectively. For UMI-counts data, we usually use 'nb'. Default is 'nb'.
#' @param nCores An integer. The number of cores to use for Parallel processing.
#' @param parallelization A string indicating the specific parallelization function to use.
#' Must be one of 'mcmapply', 'bpmapply', or 'pbmcmapply', which corresponds to the parallelization function in the package
#' \code{parallel},\code{BiocParallel}, and \code{pbmcapply} respectively. The default value is 'pbmcmapply'.
#' @param BPPARAM A \code{MulticoreParam} object or NULL. When the parameter parallelization = 'mcmapply' or 'pbmcmapply',
#' this parameter must be NULL. When the parameter parallelization = 'bpmapply',  this parameter must be one of the
#' \code{MulticoreParam} object offered by the package 'BiocParallel. The default value is NULL.
#'
#' @return The expression matrix of the synthetic null data.
#'
#' @examples
#' data(exampleCounts)
#' nullData <- constructNull(mat = exampleCounts)
#'
#' @export constructNull
constructNull <- function(mat,
                          family = "nb",
                          nCores = 1,
                          parallelization = "pbmcmapply",
                          BPPARAM = NULL) {

  sce <- SingleCellExperiment::SingleCellExperiment(list(counts = mat))
  SummarizedExperiment::colData(sce)$cell_type <- "1"
  newData <- scDesign3::scdesign3(sce,
                                  celltype = "cell_type",
                                  pseudotime = NULL,
                                  spatial = NULL,
                                  other_covariates = NULL,
                                  mu_formula = "1",
                                  sigma_formula = "1",
                                  corr_formula = "1",
                                  family_use = family,
                                  nonzerovar = FALSE,
                                  n_cores = nCores,
                                  parallelization = parallelization,
                                  important_feature = "auto",
                                  nonnegative = FALSE,
                                  copula = "gaussian",
                                  fastmvn = TRUE)
  newMat <- newData$new_count
  return(newMat)
}
