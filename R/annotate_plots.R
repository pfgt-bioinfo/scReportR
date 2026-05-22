# =============================================================================
# R/annotate_plots.R
# Annotation visualisation — all return ggplot objects.
# =============================================================================

#' UMAP coloured by SingleR labels
#'
#' @param merged    Annotated Seurat object.
#' @param ref_name  Reference name used in SingleR (e.g. "HumanPrimaryCellAtlasData").
#' @return A patchwork ggplot, or `NULL` if columns not found.
#' @export
sc_plot_singler_umap <- function(merged, ref_name) {
  prefix  <- paste0("singler_", gsub("Data$", "", ref_name), "_")
  col_lab <- paste0(prefix, "label")
  col_pru <- paste0(prefix, "pruned")
  
  if (!col_lab %in% colnames(merged@meta.data)) {
    cli::cli_alert_warning("Column not found: {col_lab}")
    return(NULL)
  }
  
  sc_plot_umap_cols(merged, cols = c(col_lab, col_pru)) |>
    sc_tabset()
}

#' SingleR score heatmap per cluster
#'
#' @param merged    Annotated Seurat object.
#' @param ref_name  Reference name.
#' @param cfg       Config list from [sc_config()].
#' @return A pheatmap object, or `NULL`.
#' @export
sc_plot_singler_heatmap <- function(merged, ref_name, cfg) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    cli::cli_alert_warning("Package 'pheatmap' not installed")
    return(NULL)
  }

  cluster_col <- paste0("clusters_res",
                        cfg$clustering$default_resolution %||% 0.6)
  prefix      <- paste0("singler_", gsub("Data$", "", ref_name), "_")
  col_lab     <- paste0(prefix, "label")
  col_score   <- paste0(prefix, "score")

  if (!col_lab %in% colnames(merged@meta.data)) return(NULL)

  score_mat <- merged@meta.data |>
    dplyr::group_by(.data[[cluster_col]], .data[[col_lab]]) |>
    dplyr::summarise(mean_score = mean(.data[[col_score]], na.rm = TRUE),
                     .groups = "drop") |>
    tidyr::pivot_wider(names_from  = dplyr::all_of(col_lab),
                       values_from = mean_score,
                       values_fill = 0) |>
    tibble::column_to_rownames(cluster_col) |>
    as.matrix()

  pheatmap::pheatmap(score_mat,
                     color        = viridis::viridis(100),
                     cluster_rows = TRUE,
                     cluster_cols = TRUE,
                     fontsize     = 8,
                     main         = paste0("SingleR scores (", ref_name, ")"))
}

#' UMAP coloured by Azimuth predictions
#'
#' @param merged  Annotated Seurat object.
#' @param level   Azimuth annotation level: `"l1"` or `"l2"`.
#' @return A patchwork ggplot.
#' @export
sc_plot_azimuth_umap <- function(merged, level = "l1") {
  col      <- paste0("predicted.celltype.", level)
  col_conf <- paste0(col, ".score")
  
  cols <- c(col, if (col_conf %in% colnames(merged@meta.data)) col_conf)
  
  if (!col %in% colnames(merged@meta.data)) {
    cli::cli_alert_warning("Column not found: {col}")
    return(NULL)
  }
  
  sc_plot_umap_cols(merged, cols = cols) |>
    sc_tabset()
}

#' UMAP coloured by manual annotation
#'
#' @param merged    Seurat object with manual annotation.
#' @param annot_col Metadata column name (default: `"celltype_manual"`).
#' @return A ggplot object.
#' @export
sc_plot_annotation <- function(merged, annot_col = "celltype_manual") {
  if (!annot_col %in% colnames(merged@meta.data)) {
    cli::cli_alert_warning("Column not found: {annot_col}")
    return(NULL)
  }
  
  sc_plot_umap_cols(merged, cols = annot_col)[[1]]
}
