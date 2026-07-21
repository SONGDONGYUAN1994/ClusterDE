#' Check Whether Simulated Null Data Outperform ClusterDE's Default Null Data
#'
#' Evaluate simulated null count matrices against ClusterDE's default null-data
#' baseline using two reference datasets, `A549_seurat` and
#' `monocyte_10x_v3_seurat`, from the ClusterDE package. The comparison uses
#' eight metrics that assess the quality of simulated null data from different
#' perspectives.
#'
#' @param new_null_data A list of simulated null count matrices generated base on
#'   the `A549_seurat` data included in the ClusterDE package. For the fairest
#'   comparison with the baseline, this should ideally contain 20 elements. One
#'   matrix is sampled at random for distributional, LISI, and silhouette
#'   checks; the full list is used to calculate null p-values in ClusterDE.
#' @param new_null_data_monocyte A list of simulated null count matrices
#'   generated base on the `monocyte_10x_v3_seurat` data included in the ClusterDE
#'   package. For the fairest comparison with the baseline, this should ideally
#'   contain 20 elements. The full list is used to calculate null p-values in
#'   ClusterDE for the monocyte differential expression check.
#' @param threshold A named numeric vector giving pass/fail thresholds for
#'   `mean`, `var`, `cor`, `lisi`, `silhouette`, `cellline`,
#'   `monocyte_marker`, and `monocyte_housekeeping`.
#' @param seed Integer seed used for reproducible sampling, PCA/UMAP, and
#'   clustering. The default is `123`.
#' @param housekeeping_gmt Local path or URL to a GMT file containing
#'   housekeeping genes. By default, the function reads the MSigDB
#'  v2026.1 human symbols GMT file. The MSigDB gene-set card for this housekeeping 
#'  set is available at <https://www.gsea-msigdb.org/gsea/msigdb/cards/HSIAO_HOUSEKEEPING_GENES>.
#' @param compare Logical. If `TRUE`, compare the calculated metrics with the
#'   values supplied in `threshold` and return the overall pass status together
#'   with a detailed comparison table. If `FALSE`, return only the calculated
#'   metrics. Defaults to `FALSE`.
#' @return A list with two elements:
#'   \describe{
#'     \item{pass}{Logical value indicating whether more than five checks pass.}
#'     \item{comparison_result}{Data frame comparing the default thresholds,
#'       observed metrics, check directions, and pass/fail results.}
#'   }
#'
#'
#' @examples
#' \dontrun{
#' library(ClusterDE)
#'
#' utils::data(A549_seurat, package = "ClusterDE")
#' utils::data(monocyte_10x_v3_seurat, package = "ClusterDE")
#'
#' RNGkind("L\'Ecuyer-CMRG")
#' set.seed(123)
#'
#' cellline_null <- constructNull(
#'   A549_seurat,
#'   nRep = 20,
#'   corrCut = 0
#' )
#'
#' monocyte_null <- constructNull(
#'   monocyte_10x_v3_seurat,
#'   nRep = 20,
#'   corrCut = 0
#' )
#'
#' threshold <- checkNullData(
#'   cellline_null,
#'   monocyte_null,
#'   compare = FALSE
#' )
#'
#' # Replace cellline_null and monocyte_null with null datasets generated
#' # from the A549 and monocyte reference datasets, respectively, using the
#' # simulator you want to evaluate.
#' df <- checkNullData(
#'   cellline_null,
#'   monocyte_null,
#'   threshold = threshold,
#'   compare = TRUE
#' )
#' }
#' @export checkNullData

checkNullData <- function(
  new_null_data,
  new_null_data_monocyte,
  threshold,
  seed = 123,
  housekeeping_gmt = "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2026.1.Hs/msigdb.v2026.1.Hs.symbols.gmt",
  compare = F
) {
  utils::data(A549_seurat, package = "ClusterDE")

  if (length(new_null_data) == 0) {
    stop("new_null_data must contain at least one null dataset.")
  }
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  idx <- sample(seq_along(new_null_data), 1)
  new_null_sample <- new_null_data[[idx]]

  ref_data <- as.matrix(A549_seurat[["RNA"]]$counts)
  new_null_sample <- as.matrix(new_null_sample)

  if (nrow(ref_data) != nrow(new_null_sample)) {
    stop("ref_data and new_null_data must have the same number of genes.")
  }
  if (!is.null(rownames(ref_data)) && !is.null(rownames(new_null_sample))) {
    if (!all(rownames(ref_data) == rownames(new_null_sample))) {
      stop("Gene names in ref_data and new_null_data are not in the same order.")
    }
  }

  ref_data_log <- log1p(ref_data)
  new_null_sample_log <- log1p(new_null_sample)

  ref_mean <- matrixStats::rowMeans2(ref_data_log)
  new_mean <- matrixStats::rowMeans2(new_null_sample_log)
  mean_cor <- stats::cor(ref_mean, new_mean, method = "pearson")

  ref_var <- matrixStats::rowVars(ref_data_log)
  new_var <- matrixStats::rowVars(new_null_sample_log)
  var_cor <- stats::cor(ref_var, new_var, method = "pearson")

  cv <- ref_var / ref_mean
  cv[!is.finite(cv)] <- NA
  top_gene_order <- order(cv, decreasing = TRUE, na.last = NA)
  top_gene_order <- top_gene_order[seq_len(min(100, length(top_gene_order)))]

  ref_gene_cor <- stats::cor(
    t(ref_data_log[top_gene_order, , drop = FALSE]),
    method = "kendall"
  )
  new_gene_cor <- stats::cor(
    t(new_null_sample_log[top_gene_order, , drop = FALSE]),
    method = "kendall"
  )
  ref_gene_cor[is.na(ref_gene_cor)] <- 0
  new_gene_cor[is.na(new_gene_cor)] <- 0

  hclust_res <- stats::hclust(stats::dist(ref_gene_cor))
  gene_order <- hclust_res$order

  ref_gene_cor <- ref_gene_cor[gene_order, gene_order]
  new_gene_cor <- new_gene_cor[gene_order, gene_order]

  ref_vec <- ref_gene_cor[upper.tri(ref_gene_cor)]
  new_vec <- new_gene_cor[upper.tri(new_gene_cor)]

  gene_gene_cor <- stats::cor(ref_vec, new_vec, method = "pearson")

  mean_lisi_res <- computeUmapLisi(ref_data, new_null_sample, seed = seed)
  mean_lisi <- mean_lisi_res$mLISI

  sil <- computeSilhouetteScore(new_null_sample, mean_lisi_res$simu_PCA, seed = seed)

  original_markers <- Seurat::FindMarkers(
    A549_seurat,
    ident.1 = 0,
    ident.2 = 1,
    min.pct = 0,
    logfc.threshold = 0
  )
  original_pval <- original_markers$p_val
  names(original_pval) <- rownames(original_markers)
  null_pval <- calcNullPval(new_null_data)
  cellline_res <- callDE(original_pval, null_pval$p, threshold = "BC")
  cellline_deg_num <- sum(cellline_res$record >= 0.5)

  utils::data(monocyte_10x_v3_seurat, package = "ClusterDE")
  original_markers <- Seurat::FindMarkers(
    monocyte_10x_v3_seurat,
    ident.1 = 0,
    ident.2 = 1,
    min.pct = 0,
    logfc.threshold = 0
  )
  original_pval <- original_markers$p_val
  names(original_pval) <- rownames(original_markers)
  null_pval <- calcNullPval(new_null_data_monocyte)
  monocyte_res <- callDE(original_pval, null_pval$p)

  utils::data(human_pbmc_marker, package = "ClusterDE")

  monocyte_markers_geneset <- union(
    setdiff(human_pbmc_marker$`CD14+ monocyte`, human_pbmc_marker$`CD16+ monocyte`),
    setdiff(human_pbmc_marker$`CD16+ monocyte`, human_pbmc_marker$`CD14+ monocyte`)
  )
  monocyte_markers_geneset <- stringr::str_remove(
    monocyte_markers_geneset,
    pattern = "\\+"
  )
  monocyte_markers_geneset <- stringr::str_remove(
    monocyte_markers_geneset,
    pattern = "\\-"
  )
  monocyte_markers_geneset <- monocyte_markers_geneset[
    monocyte_markers_geneset != "TYROBP"
  ]
  monocyte_markers_geneset <- data.frame(
    term = "CD14+/CD16+ Monocyte Markers",
    gene = monocyte_markers_geneset
  )

  if (is.null(housekeeping_gmt)) {
    stop("housekeeping_gmt must point to a local or remote GMT file.")
  }
  hkp_geneset <- gson::read.gmt(housekeeping_gmt)
  if ("HSIAO_HOUSEKEEPING_GENES" %in% hkp_geneset$term) {
    hkp_geneset <- hkp_geneset[hkp_geneset$term == "HSIAO_HOUSEKEEPING_GENES", , drop = FALSE]
  }
  if (nrow(hkp_geneset) == 0) {
    stop("No housekeeping genes were found in housekeeping_gmt.")
  }
  hkp_geneset$term <- "Housekeeping Genes"

  deg <- monocyte_res$gene[monocyte_res$record >= 0.5]

  marker_counts <- sum(deg %in% monocyte_markers_geneset$gene)
  hkp_count <- sum(deg %in% hkp_geneset$gene)

  curr_metric <- c(
    mean = round(mean_cor, 3),
    var = round(var_cor, 3),
    cor = round(gene_gene_cor, 3),
    lisi = mean_lisi,
    silhouette = sil,
    cellline = cellline_deg_num,
    monocyte_marker = marker_counts,
    monocyte_housekeeping = hkp_count
  )

  if (compare) {
    if (is.null(threshold)) {
      stop("`threshold` must be provided when `check = TRUE`.")
    }
    if (
      is.null(names(threshold)) ||
        !identical(names(threshold), names(curr_metric))
    ) {
      stop("Please provide proper values for the `threshold` parameter.")
    }
    df <- data.frame(
      ClusterDE_default = threshold,
      Current_null_data = curr_metric,
      Check = c(
        rep("Current_null_data >= ClusterDE_default", 4),
        "Current_null_data <= ClusterDE_default",
        "Current_null_data = ClusterDE_default = 0",
        "Current_null_data >= ClusterDE_default",
        "Current_null_data <= ClusterDE_default"
      ),
      Pass = c(
        curr_metric[1:4] >= threshold[1:4],
        curr_metric[5] <= threshold[5],
        curr_metric[6] == threshold[6],
        curr_metric[7] >= threshold[7],
        curr_metric[8] <= threshold[8]
      )
    )

    passed <- sum(df$Pass) > 5
    if (passed) {
      message("Null data check passed.")
    } else {
      message("Null data check failed.")
    }

    return(list(pass = passed, comparison_result = df))
  } else {
    return(curr_metric)
  }


}


#' Compute UMAP-Based LISI for Reference and Simulated Null Data
#'
#' Project reference and simulated count matrices into a shared PCA/UMAP space
#' and calculate an mean LISI score.
#'
#' @param ref_data Reference count matrix with genes in rows and cells in
#'   columns.
#' @param new_null_sample Simulated null count matrix with the same genes, in the
#'   same order, as `ref_data`.
#' @param n_pc Maximum number of principal components to use.
#' @param center Logical; passed to [stats::prcomp()].
#' @param scale. Logical; passed to [stats::prcomp()].
#' @param seed Integer seed used for reproducibility.
#'
#' @return A list with `mLISI`, the mean LISI score, and `simu_PCA`,
#'   the projected PCA coordinates for the simulated null data.
#'
#' @keywords internal
computeUmapLisi <- function(
  ref_data,
  new_null_sample,
  n_pc = 50,
  center = TRUE,
  scale. = TRUE,
  seed = 123
) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  ref_obj <- Seurat::CreateSeuratObject(counts = ref_data)
  new_obj <- Seurat::CreateSeuratObject(counts = new_null_sample)

  ref_obj <- Seurat::NormalizeData(ref_obj)
  new_obj <- Seurat::NormalizeData(new_obj)

  ref_log <- Seurat::GetAssayData(ref_obj, assay = "RNA", layer = "data")
  new_log <- Seurat::GetAssayData(new_obj, assay = "RNA", layer = "data")

  mat_ref <- t(as.matrix(ref_log))
  mat_new <- t(as.matrix(new_log))

  gene_var <- matrixStats::colVars(mat_ref)
  keep_genes <- gene_var > 0

  mat_ref <- mat_ref[, keep_genes, drop = FALSE]
  mat_new <- mat_new[, keep_genes, drop = FALSE]

  ref_pca_fit <- stats::prcomp(mat_ref, center = center, scale. = scale.)

  n_pc_use <- min(n_pc, ncol(ref_pca_fit$x))

  ref_pca <- ref_pca_fit$x[, seq_len(n_pc_use), drop = FALSE]

  new_pca <- stats::predict(
    ref_pca_fit,
    newdata = mat_new
  )[, seq_len(n_pc_use), drop = FALSE]

  ref_umap_fit <- umap::umap(ref_pca)

  ref_umap <- ref_umap_fit$layout
  colnames(ref_umap) <- c("UMAP1", "UMAP2")

  new_umap <- stats::predict(
    object = ref_umap_fit,
    data = new_pca
  )
  colnames(new_umap) <- c("UMAP1", "UMAP2")

  colnames(new_null_sample) <- paste0("Simulated_", colnames(new_null_sample))
  colnames(ref_data) <- paste0("Reference_", colnames(ref_data))

  count_combine <- cbind(new_null_sample, ref_data)

  sce_combine <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = count_combine)
  )

  SingleCellExperiment::reducedDim(sce_combine, "UMAP") <- rbind(
    new_umap,
    ref_umap
  )

  SummarizedExperiment::colData(sce_combine)$Method <- c(
    rep("Simulated", ncol(new_null_sample)),
    rep("Reference", ncol(ref_data))
  )

  sce_umap <- CellMixS::evalIntegration(
    metrics = "isi",
    sce_combine,
    k = 50,
    n_dim = 2,
    cell_min = 4,
    res_name = "weighted_isi",
    group = "Method",
    dim_red = "UMAP"
  )

  UMAP_lisi <- round(
    mean(SummarizedExperiment::colData(sce_umap)$weighted_isi, na.rm = TRUE),
    2
  )

  list(mLISI = UMAP_lisi, simu_PCA = new_pca)
}


#' Compute a Silhouette Score for Clusters in Null Data
#'
#' Cluster a simulated null count matrix with Seurat and calculate the mean
#' silhouette width from supplied PCA coordinates.
#'
#' @param count_mat Simulated null count matrix with genes in rows and cells in
#'   columns.
#' @param pca PCA coordinates for the simulated null cells.
#' @param seed Integer seed used for reproducibility.
#'
#' @return Numeric mean silhouette score rounded to two decimal places.
#'
#' @keywords internal
computeSilhouetteScore <- function(count_mat, pca, seed = 123) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)

  data <- Seurat::CreateSeuratObject(counts = count_mat)
  data <- Seurat::NormalizeData(data)
  data <- Seurat::FindVariableFeatures(data)
  data <- Seurat::ScaleData(data)
  data <- Seurat::RunPCA(data, npcs = min(50, ncol(data) - 1, nrow(data) - 1))
  data <- Seurat::FindNeighbors(data)

  right <- 0.3
  data <- Seurat::FindClusters(data, resolution = right)
  number_of_clusters <- length(unique(data$seurat_clusters))
  while (number_of_clusters < 2) {
    right <- right * 2
    data <- Seurat::FindClusters(data, resolution = right)
    number_of_clusters <- length(unique(data$seurat_clusters))
  }

  left <- 0
  while (number_of_clusters != 2) {
    mid <- (left + right) / 2
    data <- Seurat::FindClusters(data, resolution = mid)
    number_of_clusters <- length(unique(data$seurat_clusters))
    if (number_of_clusters < 2) {
      left <- mid
    } else {
      right <- mid
    }
  }

  coords <- as.matrix(pca)
  cluster_col <- paste0(Seurat::DefaultAssay(data), "_snn_res.", right)
  if (!cluster_col %in% colnames(data@meta.data)) {
    cluster_col <- "seurat_clusters"
  }
  cl <- as.numeric(data@meta.data[, cluster_col])

  sil_score <- round(
    mean(cluster::silhouette(cl, stats::dist(coords))[, "sil_width"], na.rm = TRUE),
    2
  )

  sil_score
}
