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
#' @param fastVersion A logic value. If TRUE, the fast approximation is used. Default is FALSE.
#' @param corrCut A numeric value. The cutoff for non-zero proportions in genes used in modelling correlation.
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
                          fastVersion = FALSE,
                          corrCut = 0.2,
                          BPPARAM = NULL) {
  isSparse <- methods::is(mat, "sparseMatrix")

  if(!fastVersion) {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = mat))
    SummarizedExperiment::colData(sce)$cell_type <- "1"
    newData <- scDesign3::scdesign3(sce,
                                    celltype = "cell_type",
                                    pseudotime = NULL,
                                    spatial = NULL,
                                    other_covariates = NULL,
                                    empirical_quantile = FALSE,
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
                                    fastmvn = FALSE)
    newMat <- newData$new_count
  } else {
    tol <- 1e-5
    mat <- as.matrix(mat)
    n_gene <- dim(mat)[1]
    n_cell <- dim(mat)[2]
    gene_names <- rownames(mat)

    qc <- apply(mat, 1, function(x){
      return(length(which(x < tol)) > length(x) - 3)
    })
    if(length(which(qc)) == 0){
      filtered_gene <- NULL
    }else{
      filtered_gene <- names(which(qc))
    }

    mat_filtered <- mat[!qc, ]

    important_feature <- names(which(matrixStats::rowMeans2(mat_filtered!=0) > corrCut))

    unimportant_feature <- setdiff(gene_names, union(important_feature, filtered_gene))

    para_feature <- rownames(mat_filtered)

    mat_corr <- t(mat_filtered[important_feature, ])

    if(family == "nb") {
      para <- parallel::mclapply(X = seq_len(dim(mat_filtered)[1]),
                                 FUN = function(x) {
                                   tryCatch({
                                     res <- fitdistrplus::fitdist(mat_filtered[x, ], "nbinom", method = "mle")$estimate
                                     res},
                                     error = function(cond) {
                                       message(paste0(x, "is problematic with NB MLE; using Poisson MME instead."))
                                       fit_para <- fitdistrplus::fitdist(mat_filtered[x, ], "pois", method = "mme")$estimate
                                       res <- c(NA, fit_para)
                                       names(res) <- c("size", "mu")
                                       res
                                     })
                                 },
                                 mc.cores = nCores)
      para <- t(simplify2array(para))
      rownames(para) <- para_feature
    } else if (family == "poisson") {
      para <- parallel::mclapply(X = seq_len(dim(mat_filtered)[1]),
                                 FUN = function(x) {
                                   tryCatch({
                                     res <- fitdistrplus::fitdist(mat_filtered[x, ], "pois", method = "mle")$estimate
                                     res},
                                     error = function(cond) {
                                       message(paste0(x, "is problematic with Poisson MLE; using Poisson MME instead."))
                                       fit_para <- fitdistrplus::fitdist(mat_filtered[x, ], "pois", method = "mme")$estimate
                                       #res <- c(NA, mu = fit_para)
                                       #names(res) <- c("size", "mu")
                                       res
                                     })
                                 },
                                 mc.cores = nCores)
      para <- simplify2array(para)
      names(para) <- para_feature
    } else {
      stop("FastVersion only supports NB and Poisson.")
    }


    p_obs <- rvinecopulib::pseudo_obs(mat_corr)
    normal_obs <- stats::qnorm(p_obs)

    corrlation <- function(x) {
      mat <- t(x) - matrixStats::colMeans2(x)
      mat <- mat / sqrt(matrixStats::rowSums2(mat^2))
      tcrossprod(mat)
    }

    corr_mat <- corrlation(normal_obs)
    diag(corr_mat) <- diag(corr_mat) + tol
    new_mvn <- mvnfast::rmvn(n_cell,
                  mu = rep(0, dim(corr_mat)[1]),
                  sigma = corr_mat,
                  isChol = FALSE,
                  ncores = nCores)
    colnames(new_mvn) <- important_feature
    new_mvp <- stats::pnorm(new_mvn)

    newMat <- matrix(0, nrow = n_gene, ncol = n_cell)
    rownames(newMat) <- gene_names
    colnames(newMat) <- paste0("Cell", seq_len(n_cell))

    if(length(unimportant_feature) > 0) {
        unimportant_mat <- parallel::mclapply(unimportant_feature, function(x) {
          if(family == "nb") {
            if(is.na(para[x, 1])) {
              stats::rpois(n = n_cell, lambda = para[x, 2])
            } else {
              stats::rnbinom(n = n_cell, size = para[x, 1], mu = para[x, 2])
            }
            stats::rnbinom(n = n_cell, size = para[x, 1], mu = para[x, 2])
          } else {
            stats::rpois(n = n_cell, lambda = para[x])
          }
        }, mc.cores = nCores)

        unimportant_mat <- t(simplify2array(unimportant_mat))
        rownames(unimportant_mat) <- unimportant_feature

        newMat[unimportant_feature, ] <- unimportant_mat
    }

    important_mat <- parallel::mclapply(important_feature, function(x) {
      if(family == "nb") {
        stats::qnbinom(p = as.vector(new_mvp[, x]), size = para[x, 1], mu = para[x, 2])}
      else {
        stats::qpois(p = as.vector(new_mvp[, x]), lambda = para[x])
      }
    }, mc.cores = nCores)

    important_mat <- t(simplify2array(important_mat))
    rownames(important_mat) <- important_feature

    newMat[important_feature, ] <- important_mat
    if(isSparse){
      newMat <- Matrix::Matrix(newMat, sparse = TRUE)
    }
  }
  return(newMat)
}