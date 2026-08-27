#' Construct the synthetic null data
#'
#' \code{constructNull} takes the target data as the input and returns the corresponding synthetic null data.
#'
#' This function constructs the synthetic null data based on the target data (real data). The input is a expression matrix (gene by cell); the user should specify a distribution, which is usually Negative Binomial for count matrix.
#'
#' @param obj A Seurat object. The reference data.
#' @param approximation A string of either "none", "fast" or "pca". For a high-latitude scRNA counting matrix, use "fast" as approximation can increase the speed of data generation while ensuring accuracy. For high-dimensional scRNA data (gene number is much larger than cell number), use "pca" as approximation. Default is "none".
#' @param corr_cut A numeric value. The cutoff for non-zero proportions in genes used in modelling correlation. Default is 0.1. All features will be used if approximation is set to "pca", or when data_type is "bulk_microarray".
#' @param family A string of the distribution of your data.
#' Must be one of 'nb', 'binomial', 'poisson', 'zip', 'zinb' or 'gaussian', which represent 'poisson distribution',
#' 'negative binomial distribution', 'zero-inflated poisson distribution', 'zero-inflated negative binomail distribution',
#' and 'gaussian distribution' respectively. For UMI-counts data, we usually use 'nb'. For bulk microarray data, use 'gaussian'. Default is 'nb'.
#' @param data_type A string of either "scRNA", "scATAC", "spatial", "cellline", "microbiome" or "bulk_microarray". Default is "scRNA".
#' @param formula A string of the mu parameter formula. It defines the relationship between gene expression in synthetic null data and the extra covariates. Default is 1 (cell type case).
#' For example, if your input data is a spatial data with X, Y coordinates, the formula can be 's(X, Y, bs = 'gp', k = 4)'.
#' @param if_sparse A logic value. For high-dimensional data (gene number is much larger than cell number), if a sparse correlation estimation will be used. Default is FALSE.
#' @param n_cores An integer. The number of cores to use for Parallel processing.
#' @param n_pcs A numeric value. Number of PCs to use when usePca=T. Default is 200.
#' @param n_rep An integer. The number of sampled synthetic null datasets. Default value is 1.
#' @param other_covariates A list of the extra covariates used in \code{formula}. For example, the 2D spatial coordinates. Default is NULL.
#' @param seed Random seed. Default is 123
#'
#' @return The expression matrix of the synthetic null data.
#'
#' @importFrom gamlss.dist dZIP pZIP qZIP rZIP ZIP
#' @export constructNull
constructNull <- function(
  obj,
  approximation = "none",
  corr_cut = 0.1,
  data_type = "scRNA",
  family = "nb",
  formula = "1",
  if_sparse = F,
  n_cores = 1,
  n_pcs = 200,
  n_rep = 1,
  other_covariates = NULL,
  seed = 123
) {
  mat <- Seurat::GetAssayData(obj, layer = "counts")
  ## Check if we should use sparse matrix.
  isSparse <- methods::is(mat, "sparseMatrix")

  if (n_rep < n_cores) {
    n_cores <- n_rep
  }

  supported_data_types <- c("scRNA", "scATAC", "spatial", "cellline", "microbiome", "bulk_microarray")
  if (!data_type %in% supported_data_types) {
    stop(sprintf(
      "Invalid data_type: '%s'. Supported types are: %s",
      data_type,
      paste(supported_data_types, collapse = ", ")
    ))
  }

  # Add fake variable for cell type constraint
  obj@meta.data$fake_variable <- 1

  synthetic_null_list <- if (approximation == "pca") {
    # Construct PCA
    message("Contruct PCA")
    non_zero_genes <- apply(mat, 1, var) != 0
    n_zero_genes <- sum(!non_zero_genes)
    if (n_zero_genes > 0) {
      message(n_zero_genes, " genes removed due to zero variances")
    }
    mat <- mat[non_zero_genes, , drop = F]
    normalized_mat <- t(logcp10k(as.matrix(mat)))
    pca_res <- prcomp(
      normalized_mat,
      center = T,
      scale. = T
    )
    pca_loading <- pca_res$rotation
    rownames(pca_loading) <- rownames(mat)
    pca_score <- pca_res$x
    rownames(pca_score) <- colnames(mat)
    ## get the bootstrapped residuals
    reconstructed_mat <- pca_score[, 1:n_pcs] %*% t(pca_loading[, 1:n_pcs])
    pca_intput <- sweep(normalized_mat, 2, pca_res$center, "-")
    pca_intput <- sweep(pca_intput, 2, pca_res$scale, "/")
    residuals <- pca_intput - reconstructed_mat

    pca_sce <- SingleCellExperiment::SingleCellExperiment(list(counts = t(pca_score[, 1:n_pcs])), colData = obj@meta.data)

    set.seed(seed)
    message("Construct scDesign3 data")
    data <- scDesign3::construct_data(
      sce = pca_sce,
      assay_use = "counts",
      celltype = "fake_variable",
      pseudotime = NULL,
      spatial = other_covariates,
      other_covariates = NULL,
      corr_by = "ind"
    )

    message("Fit marginal")
    marginal <- scDesign3::fit_marginal(
      data = data,
      predictor = "gene",
      mu_formula = formula,
      sigma_formula = "1",
      family_use = "gaussian",
      n_cores = n_cores,
      parallelization = "mapply"
    )

    message("Fit copula")
    copula <- scDesign3::fit_copula(
      sce = pca_sce,
      assay_use = "counts",
      input_data = data$dat,
      marginal_list = marginal,
      family_use = "gaussian",
      n_cores = n_cores,
      parallelization = "mapply"
    )

    para_list <- scDesign3::extract_para(
      sce = pca_sce,
      assay_use = "counts",
      marginal_list = marginal,
      family_use = "gaussian",
      new_covariate = data$newCovariate,
      data = data$dat,
      n_cores = n_cores,
      parallelization = "mapply"
    )

    message(paste0("Generate null data of ", n_rep, " replicates"))
    new_count_list <- suppressMessages(bettermc::mclapply(1:n_rep, function(b) {
      set.seed(seed + b)
      new_count <- scDesign3::simu_new(
        sce = pca_sce,
        assay_use = "counts",
        mean_mat = para_list$mean_mat,
        sigma_mat = para_list$sigma_mat,
        zero_mat = para_list$zero_mat,
        quantile_mat = NULL,
        copula_list = copula$copula_list,
        n_cores = 1,
        family_use = "gaussian",
        nonnegative = FALSE,
        nonzerovar = FALSE,
        input_data = data$dat,
        new_covariate = data$newCovariate,
        important_feature = copula$important_feature,
        filtered_gene = data$filtered_gene
      )

      new_mat <- t(new_count) %*% t(pca_loading[, 1:n_pcs])
      residuals_bootstrap <- apply(residuals, 2, function(x) sample(x, length(x), replace = TRUE))
      rownames(residuals_bootstrap) <- rownames(residuals)
      new_mat <- new_mat + residuals_bootstrap
      new_mat <- sweep(new_mat, 2, pca_res$scale, `*`)
      new_mat <- sweep(new_mat, 2, pca_res$center, `+`)
      t(new_mat)
    }, mc.cores = n_cores, mc.retry = 5))
    new_count_list
  } else if (
    (data_type == "scRNA") ||
      (data_type == "scATAC") ||
      (data_type == "cellline")
  ) {
    set.seed(seed)
    tol <- 1e-5
    mat <- as.matrix(mat)
    n_gene <- dim(mat)[1]
    n_cell <- dim(mat)[2]
    gene_names <- rownames(mat)

    qc <- apply(mat, 1, function(x) {
      return(length(which(x < tol)) > length(x) - 3)
    })
    if (length(which(qc)) == 0) {
      filtered_gene <- NULL
    } else {
      filtered_gene <- names(which(qc))
      message(
        paste0(
          length(which(qc)),
          " genes have no more than 2 non-zero values; ignore fitting and return all 0s."
        )
      )
    }

    mat_filtered <- mat[!qc,]
    para_feature <- rownames(mat_filtered)

    ## Marginal fitting

    if (family == "nb") {
      para <- parallel::mclapply(
        X = seq_len(dim(mat_filtered)[1]),
        FUN = function(x) {
          tryCatch({
            res <- suppressWarnings(fitdistrplus::fitdist(mat_filtered[x,], "nbinom", method = "mle")$estimate)
            res
          }, error = function(cond) {
            message(paste0(
              x,
              " is problematic with NB MLE; using Poisson MME instead."
            ))
            fit_para <- suppressWarnings(fitdistrplus::fitdist(mat_filtered[x,], "pois", method = "mme")$estimate)
            res <- c(NA, fit_para)
            names(res) <- c("size", "mu")
            res
          })
        },
        mc.cores = n_cores
      )
      para <- t(simplify2array(para))
      rownames(para) <- para_feature

      if (sum(is.na(para[, 2])) > 0) {
        warning("NA produces in mean estimate; using 0 instead.")
        para[, 2][is.na(para[, 2])] <- 0
      }

    }
    else if (family == "poisson") {
      para <- parallel::mclapply(
        X = seq_len(dim(mat_filtered)[1]),
        FUN = function(x) {
          tryCatch({
            res <- fitdistrplus::fitdist(mat_filtered[x,], "pois", method = "mle")$estimate
            res
          }, error = function(cond) {
            message(paste0(
              x,
              "is problematic with Poisson MLE; using Poisson MME instead."
            ))
            res
          })
        },
        mc.cores = n_cores
      )
      para <- simplify2array(para)
      names(para) <- para_feature
      if (sum(is.na(para)) > 0) {
        warning("NA produces in mean estimate; using 0 instead.")
        para[is.na(para)] <- 0
      }
    } else if (family == "zip") {
      para <- parallel::mclapply(
        X = seq_len(dim(mat_filtered)[1]),
        FUN = function(x) {
          tryCatch({
            res <- suppressWarnings(
              fitdistrplus::fitdist(
                mat_filtered[x,],
                "ZIP",
                method = "mle",
                start = list(mu = mean(mat[x,]), sigma = 0.1)
              )$estimate
            )
            res
          }, error = function(cond) {
            message(paste0(
              x,
              " is problematic with NB MLE; using Poisson MME instead."
            ))
            fit_para <- suppressWarnings(fitdistrplus::fitdist(mat_filtered[x,], "pois", method = "mme")$estimate)
            res <- c(fit_para, NA)
            names(res) <- c("mu", "sigma")
            res
          })
        },
        mc.cores = n_cores
      )
      para <- t(simplify2array(para))
      rownames(para) <- para_feature

      if (sum(is.na(para[, 1])) > 0) {
        warning("NA produces in mean estimate; using 0 instead.")
        para[, 1][is.na(para[, 1])] <- 0
      }
    } else {
      stop("scRNA data distribution family should be one of 'nb', 'poisson' or 'zip'.")
    }

    ## Now we get the para matrix. You can modify it here. First column is the dispersion and second column is the mean.

    ## Copula fitting
    important_feature <- names(which(rowMeans(mat_filtered != 0) > corr_cut))

    if (length(important_feature) > 1) {
      unimportant_feature <- setdiff(gene_names, union(important_feature, filtered_gene))

      mat_corr <- t(mat_filtered[important_feature,])
      corr_prop <- round(length(important_feature) / n_gene, 3)
      p_obs <- rvinecopulib::pseudo_obs(mat_corr)
      normal_obs <- stats::qnorm(p_obs)

      message(paste0(
        corr_prop * 100,
        "% of genes are used in correlation modelling."
      ))

      if (if_sparse) {
        corr_mat <- scDesign3::sparse_cov(
          normal_obs,
          method = 'qiu',
          operator = 'hard',
          corr = TRUE
        )
      } else {
        corr_mat <- coop::pcor(normal_obs)
      }

      diag(corr_mat) <- diag(corr_mat) + tol

      ####
      if (approximation == "none") {
        #get parameters for Cholesky decomposition factor
        cdf <- chol(corr_mat)
      } else {
        # get parameters for block sampling

        #It is guaranteed that non-positive definite matrices can also be Cholesky decomposed
        approx_chol_eigen_direct <- function(mat,
                                             eps = 1e-6,
                                             verbose = TRUE) {
          e <- eigen(mat, symmetric = TRUE)
          if (any(e$values < eps)) {
            if (verbose)
              message("Eigenvalue correction applied.")
            e$values[e$values < eps] <- eps
          }
          sqrt_vals <- sqrt(e$values)
          chol_factor <- e$vectors %*% diag(sqrt_vals)
          return(chol_factor)
        }

        simple_block_chol <- function(mat, eps = 1e-6) {
          #Positive definiteness correction for processing matrices in blocks
          k <- nrow(mat)
          idx <- split(1:k, cut(1:k, breaks = 4, labels = FALSE))

          L_blocks <- lapply(idx, \(i) approx_chol_eigen_direct(mat[i, i], eps))
          L <- as.matrix(Matrix::bdiag(L_blocks))

          diag(L) <- diag(L) * (1 + eps)
          L
        }

        d <- nrow(corr_mat)
        k <- ceiling(d / 2)

        L12 <- corr_mat[1:k, (k + 1):d]
        svd_res <- svd(L12)
        d_all <- svd_res$d
        r <- sum(d_all > 1e-3)

        U <- svd_res$u[, 1:r]
        V <- svd_res$v[, 1:r]
        D_root <- sqrt(d_all[1:r])

        U_t <- t(U * outer(rep(1, nrow(U)), D_root))

        V_t <- t(V * outer(rep(1, nrow(V)), D_root))

        L11 <- corr_mat[1:k, 1:k] - crossprod(U_t)
        L22 <- corr_mat[(k + 1):d, (k + 1):d] - crossprod(V_t)

        l_bm11 <- simple_block_chol(L11) #ensure positive definition
        l_bm22 <- simple_block_chol(L22) #ensure positive definition


        block_mvn_sample <- function(n_cell,
                                     l_bm11,
                                     l_bm22,
                                     U_t,
                                     V_t,
                                     k,
                                     d,
                                     ncores = ncores) {
          X <- matrix(0, nrow = n_cell, ncol = d)
          X[, 1:k] <- mvnfast::rmvn(
            n_cell,
            mu = rep(0, k),
            sigma = l_bm11,
            isChol = TRUE,
            ncores = ncores
          )
          X[, (k + 1):d] <- mvnfast::rmvn(
            n_cell,
            mu = rep(0, d - k),
            sigma = l_bm22,
            isChol = TRUE,
            ncores = ncores
          )


          if (!is.null(U_t)) {
            r <- nrow(U_t)
            Z <- matrix(zigg::zrnorm(n_cell * r), nrow = n_cell)

            X[, 1:k] <- X[, 1:k] + matrix_multiplication_cpp(Z, U_t)
            X[, (k + 1):d] <- X[, (k + 1):d] + matrix_multiplication_cpp(Z, V_t)

          }
          return(X)
        }

      }
      ## Start sampling
      if (n_rep == 1) {
        if (approximation == "none") {
          new_mvn <- mvnfast::rmvn(
            n_cell,
            mu = rep(0, dim(corr_mat)[1]),
            sigma = cdf,
            isChol = TRUE,
            ncores = n_cores
          )
        } else {
          new_mvn <- block_mvn_sample(
            n_cell = n_cell,
            l_bm11 = l_bm11,
            l_bm22 = l_bm22,
            U_t = U_t,
            V_t = V_t,
            k = k,
            d = d,
            ncores = n_cores
          )

        }

        colnames(new_mvn) <- important_feature
        new_mvp <- stats::pnorm(new_mvn)

        newMat <- matrix(0, nrow = n_gene, ncol = n_cell)
        rownames(newMat) <- gene_names
        colnames(newMat) <- paste0("Cell", seq_len(n_cell))

        if (length(unimportant_feature) > 0) {
          unimportant_mat <- parallel::mclapply(unimportant_feature, function(x) {
            if (family == "nb") {
              if (is.na(para[x, 1])) {
                stats::rpois(n = n_cell, lambda = para[x, 2])
              } else {
                stats::rnbinom(n = n_cell,
                               size = para[x, 1],
                               mu = para[x, 2])
              }
            } else if (family == "poisson") {
              stats::rpois(n = n_cell, lambda = para[x])
            } else if (family == "zip") {
              if (is.na(para[x, 2])) {
                stats::rpois(n = n_cell, lambda = para[x, 1])
              } else {
                rZIP(n = n_cell,
                     sigma = para[x, 2],
                     mu = para[x, 1])
              }
            } else {
              stop("Family must be in nb, poisson, or zip.")
            }
          }, mc.cores = n_cores)

          unimportant_mat <- t(simplify2array(unimportant_mat))
          rownames(unimportant_mat) <- unimportant_feature

          newMat[unimportant_feature,] <- unimportant_mat
        }

        important_mat <- parallel::mclapply(important_feature, function(x) {
          if (family == "nb") {
            if (is.na(para[x, 1])) {
              stats::qpois(p = as.vector(new_mvp[, x]), lambda = para[x, 2])
            } else {
              stats::qnbinom(p = as.vector(new_mvp[, x]),
                             size = para[x, 1],
                             mu = para[x, 2])
            }
          } else if (family == "poisson") {
            stats::qpois(p = as.vector(new_mvp[, x]), lambda = para[x])
          } else if (family == "zip") {
            if (is.na(para[x, 2])) {
              stats::qpois(p = as.vector(new_mvp[, x]), lambda = para[x, 1])
            } else {
              qZIP(
                p = as.vector(new_mvp[, x]),
                sigma = para[x, 2],
                mu = para[x, 1]
              )
            }
          } else {
            stop("Family must be in nb, poisson, or zip.")
          }
        }, mc.cores = n_cores)

        important_mat <- t(simplify2array(important_mat))
        rownames(important_mat) <- important_feature

        newMat[important_feature,] <- important_mat
        newMat[is.na(newMat)] <- 0
        if (isSparse) {
          newMat <- Matrix::Matrix(newMat, sparse = TRUE)
        }
        newMat

      } else {
        parallel::mclapply(seq_len(n_rep), function(x) {
          if (approximation == "none") {
            new_mvn <- mvnfast::rmvn(
              n_cell,
              mu = rep(0, dim(corr_mat)[1]),
              sigma = cdf,
              isChol = TRUE,
              ncores = 1
            )
          } else {
            new_mvn <- block_mvn_sample(
              n_cell = n_cell,
              l_bm11 = l_bm11,
              l_bm22 = l_bm22,
              U_t = U_t,
              V_t = V_t,
              k = k,
              d = d,
              ncores = 1
            )

          }

          colnames(new_mvn) <- important_feature
          new_mvp <- stats::pnorm(new_mvn)

          newMat <- matrix(0, nrow = n_gene, ncol = n_cell)
          rownames(newMat) <- gene_names
          colnames(newMat) <- paste0("Cell", seq_len(n_cell))

          if (length(unimportant_feature) > 0) {
            unimportant_mat <- lapply(unimportant_feature, function(x) {
              if (family == "nb") {
                if (is.na(para[x, 1])) {
                  stats::rpois(n = n_cell, lambda = para[x, 2])
                } else {
                  stats::rnbinom(n = n_cell,
                                 size = para[x, 1],
                                 mu = para[x, 2])
                }
              } else if (family == "poisson") {
                stats::rpois(n = n_cell, lambda = para[x])
              } else if (family == "zip") {
                if (is.na(para[x, 2])) {
                  stats::rpois(n = n_cell, lambda = para[x, 1])
                } else {
                  rZIP(n = n_cell,
                       sigma = para[x, 2],
                       mu = para[x, 1])
                }
              } else {
                stop("Family must be in nb, poisson, or zip.")
              }
            })

            unimportant_mat <- t(simplify2array(unimportant_mat))
            rownames(unimportant_mat) <- unimportant_feature

            newMat[unimportant_feature,] <- unimportant_mat
          }

          important_mat <- lapply(important_feature, function(x) {
            if (family == "nb") {
              if (is.na(para[x, 1])) {
                stats::qpois(p = as.vector(new_mvp[, x]), lambda = para[x, 2])
              } else {
                stats::qnbinom(p = as.vector(new_mvp[, x]),
                               size = para[x, 1],
                               mu = para[x, 2])
              }
            } else if (family == "poisson") {
              stats::qpois(p = as.vector(new_mvp[, x]), lambda = para[x])
            } else if (family == "zip") {
              if (is.na(para[x, 2])) {
                stats::qpois(p = as.vector(new_mvp[, x]), lambda = para[x, 1])
              } else {
                qZIP(
                  p = as.vector(new_mvp[, x]),
                  sigma = para[x, 2],
                  mu = para[x, 1]
                )
              }
            } else {
              stop("Family must be in nb, poisson, or zip.")
            }
          })

          important_mat <- t(simplify2array(important_mat))
          rownames(important_mat) <- important_feature

          newMat[important_feature,] <- important_mat
          newMat[is.na(newMat)] <- 0
          if (isSparse) {
            newMat <- Matrix::Matrix(newMat, sparse = TRUE)
          }
          newMat
        }, mc.cores = n_cores)
      }

    } else {
      message("No correlation structure. All features are independent.")
      lapply(seq_len(n_rep), function(x) {
        newMat <- matrix(0, nrow = n_gene, ncol = n_cell)
        rownames(newMat) <- gene_names
        colnames(newMat) <- paste0("Cell", seq_len(n_cell))

        para_mat <- parallel::mclapply(para_feature, function(x) {
          if (family == "nb") {
            if (is.na(para[x, 1])) {
              stats::rpois(n = n_cell, lambda = para[x, 2])
            } else {
              stats::rnbinom(n = n_cell,
                             size = para[x, 1],
                             mu = para[x, 2])
            }
            stats::rnbinom(n = n_cell,
                           size = para[x, 1],
                           mu = para[x, 2])
          } else if (family == "poisson") {
            stats::rpois(n = n_cell, lambda = para[x])
          } else if (family == "zip") {
            if (is.na(para[x, 2])) {
              stats::rpois(n = n_cell, lambda = para[x, 1])
            } else {
              gamlss.dist::rZIP(n = n_cell,
                                sigma = para[x, 2],
                                mu = para[x, 1])
            }
          } else {
            stop("Family must be in nb, poisson, or zip.")
          }
        }, mc.cores = n_cores)

        para_mat <- t(simplify2array(para_mat))
        newMat[para_feature,] <- para_mat
        newMat[is.na(newMat)] <- 0
        if (isSparse) {
          newMat <- Matrix::Matrix(newMat, sparse = TRUE)
        }
        newMat
      })
    }
  } else {
    sce <- SingleCellExperiment::SingleCellExperiment(list(counts = mat))
    SummarizedExperiment::colData(sce) <- S4Vectors::DataFrame(obj@meta.data)

    set.seed(seed)
    newData <- scDesign3::scdesign3(
      sce,
      celltype = "fake_variable",
      pseudotime = NULL,
      spatial = if (data_type == "spatial") other_covariates else NULL,
      other_covariates = if (data_type != "spatial") other_covariates else NULL,
      empirical_quantile = FALSE,
      mu_formula = formula,
      sigma_formula = "1",
      corr_formula = "1",
      family_use = family,
      nonzerovar = FALSE,
      n_cores = n_cores,
      parallelization = "mcmapply",
      important_feature = if (data_type != "bulk_microarray") corr_cut else "all",
      nonnegative = FALSE,
      copula = "gaussian",
      if_sparse = if_sparse,
      fastmvn = FALSE,
      n_rep = n_rep
    )
    newData$new_count
  }

  if (length(synthetic_null_list) == 1) {
    synthetic_null_list[[1]]
  } else {
    synthetic_null_list
  }
}
