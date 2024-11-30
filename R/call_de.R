#' Call differentially expressed (DE) genes by Clipper algorithm
#'
#' \code{callDE} takes two vectors representing each gene's significance. They are usually p-values from the target data and synthetic null data.
#'
#' This function constructs the contrast scores by taking the difference between target DE scores and null DE scores.
#'
#' @param targetScores A named numeric vector of the DE scores from the target data, e.g., the p-values between two clusters from the real data.
#' @param nullScores A named numeric vector of the DE scores from the synthetic null data, e.g., the p-values between two clusters from the null data.
#' @param nlogTrans A logical value. If the input scores are p-values, take the \code{-log10} transformation since Clipper require larger scores represent more significant DE. Default is TRUE.
#' @param FDR A numeric value of the target False Discovery Rate (FDR). Must be 'diff' or 'max'.
#' @param contrastScore A string value of the way to construct contrast scores. The choice can be
#' @param correct A logical value. If TRUE, perform the correction to make the distribution of contrast scores approximately symmetric. Default is FALSE.
#' @param threshold A string value of the threshold method. Must be 'BC' or 'DS'.
#' @param ordering A logic value. If TRUE, order the genes in the returned table by their significance. Default is TRUE.
#'
#' @return A list of target FDR, DE genes, and the detailed summary table.
#'
#' @examples
#' targetScores <- runif(10000)
#' nullScores <- runif(10000)
#' names(targetScores) <- names(nullScores) <- paste0("Gene", 1:10000)
#' res <- callDE(targetScores, nullScores, correct = FALSE)
#'
#' @export callDE
callDE <- function(targetScores,
                   nullScores,
                   nlogTrans = TRUE,
                   FDR = 0.05,
                   contrastScore = "diff",
                   correct = FALSE,
                   threshold = "BC",
                   ordering = TRUE) {
  value <- null <- target <- cs <- NULL

  if(is.null(names(targetScores))|is.null(names(nullScores))) {
    stop("Both scores should have gene names!")
  }

  tbl_target <- dplyr::as_tibble(targetScores, rownames = "Gene")
  tbl_target <- dplyr::rename(tbl_target, target = value)
  tbl_null <- dplyr::as_tibble(nullScores, rownames = "Gene")
  tbl_null <- dplyr::rename(tbl_null, null = value)

  tbl_merge <- dplyr::left_join(tbl_target, tbl_null, by = "Gene")
  if(nlogTrans) {
    tbl_merge <- dplyr::mutate(tbl_merge, target = -log10(target), null = -log10(null))
  }

  if(contrastScore == 'diff') {
    tbl_merge <- dplyr::mutate(tbl_merge, cs = target - null) ## Diff contrast scores

    if(correct) {
      if(PairedData::yuen.t.test(x = tbl_merge$target, y = tbl_merge$null, alternative = "greater", paired = TRUE, tr = 0.1)$p.value < 0.001) {
        fit <- MASS::rlm(tbl_merge$target ~ tbl_merge$null, maxit = 100)
        tbl_merge$cs <- fit$residuals
      }
    }
  } else if (contrastScore == 'max') {
    tbl_merge <- dplyr::mutate(tbl_merge, cs = max(target, null)*sign(target - null))
    } else stop("Contrast score must be constructed by 'diff' or 'max' method.")


  tbl_merge <- dplyr::mutate(tbl_merge, q = cs2q(contrastScore = tbl_merge$cs, threshold = threshold))
  if(ordering) {
    tbl_merge <- dplyr::arrange(tbl_merge, dplyr::desc(cs))
  }

  DEgenes <- as.vector(dplyr::filter(tbl_merge, q <= FDR)$Gene)
  return(list(targetFDR = FDR, DEgenes = DEgenes, summaryTable = tbl_merge))
}

cs2q <- function(contrastScore, nnull = 1, threshold = "BC"){
  stopifnot(threshold == "BC" | threshold == "DS")

  contrastScore[is.na(contrastScore)] = 0 # impute missing contrast scores with 0
  c_abs = abs(contrastScore[contrastScore != 0])
  c_abs  = sort(unique(c_abs))

  i = 1
  emp_fdp = rep(NA, length(c_abs))
  emp_fdp[1] = 1

  if(threshold == "BC") {
    while(i <= length(c_abs)){
      # print(i)
      t = c_abs[i]
      emp_fdp[i] = min((1/nnull + 1/nnull * sum(contrastScore <= -t))/ sum(contrastScore >= t),1)

      if (i >=2){emp_fdp[i] = min(emp_fdp[i], emp_fdp[i-1])}
      i = i + 1
    }
  } else if (threshold == "DS") {
    while(i <= length(c_abs)){
      # print(i)
      t = c_abs[i]
      emp_fdp[i] = min((1/nnull * sum(contrastScore <= -t))/ sum(contrastScore >= t),1)

      if (i >=2){emp_fdp[i] = min(emp_fdp[i], emp_fdp[i-1])}
      i = i + 1
    }
  } else {stop("Method must be BC or DS!")}



  c_abs = c_abs[!is.na(emp_fdp)]
  emp_fdp = emp_fdp[!is.na(emp_fdp)]
  q <- emp_fdp[match(contrastScore, c_abs)]
  q[which(is.na(q))] = 1

  return(q)
}

