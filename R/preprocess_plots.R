# =============================================================================
# R/preprocess_plots.R
# Preprocessing visualisation — all return ggplot objects.
# =============================================================================

#' Bar chart: % doublets per sample
#'
#' @param obj_list  Named list of Seurat objects after [sc_doublets()].
#' @param cfg       Config list from [sc_config()].
#' @param fill_by   Metadata column for fill colour. Defaults to best available.
#' @return A ggplot object.
#' @export
sc_plot_doublets <- function(obj_list, cfg, fill_by = NULL) {

  meta <- dplyr::bind_rows(
    lapply(obj_list, function(o) as.data.frame(o@meta.data))
  )
  fill_col    <- fill_by %||% .best_meta_col(meta)
  fill_colors <- .resolve_palette(fill_col, cfg)

  p <- meta |>
    dplyr::count(sample_id, .data[[fill_col]], scDblFinder_class) |>
    dplyr::group_by(sample_id) |>
    dplyr::mutate(pct = 100 * n / sum(n)) |>
    dplyr::filter(scDblFinder_class == "doublet") |>
    ggplot2::ggplot(ggplot2::aes(x = sample_id, y = pct,
                                 fill = .data[[fill_col]])) +
    ggplot2::geom_col(alpha = 0.85) +
    ggplot2::labs(title = "% doublets detected per sample",
                  x = NULL, y = "% doublets") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x    = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "bottom")

  if (!is.null(fill_colors))
    p <- p + ggplot2::scale_fill_manual(values = fill_colors)
  p
}

#' Violin: decontX contamination scores per sample
#'
#' @param obj_list   Named list of Seurat objects after [sc_decontX()].
#' @param cfg        Config list from [sc_config()].
#' @param threshold  Contamination score threshold to draw as a reference line.
#' @param fill_by    Metadata column for fill colour. Defaults to best available.
#' @return A ggplot object.
#' @export
sc_plot_decontX <- function(obj_list, cfg,
                             threshold = cfg$decontX$contamination_threshold,
                             fill_by   = NULL) {

  meta <- dplyr::bind_rows(
    lapply(obj_list, function(o) as.data.frame(o@meta.data))
  ) |>
    dplyr::filter(!is.na(decontX_contamination))

  if (nrow(meta) == 0L) {
    return(
      ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "decontX not run") +
        ggplot2::theme_void()
    )
  }

  fill_col    <- fill_by %||% .best_meta_col(meta)
  fill_colors <- .resolve_palette(fill_col, cfg)

  p <- ggplot2::ggplot(meta,
                       ggplot2::aes(x    = sample_id,
                                    y    = decontX_contamination,
                                    fill = .data[[fill_col]])) +
    ggplot2::geom_violin(scale = "width", alpha = 0.8) +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed",
                        color = "red") +
    ggplot2::labs(title    = "decontX contamination score",
                  subtitle = paste0("Red line = threshold (", threshold, ")"),
                  x = NULL, y = "Contamination score") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x    = ggplot2::element_text(angle = 45, hjust = 1),
                   legend.position = "bottom")

  if (!is.null(fill_colors))
    p <- p + ggplot2::scale_fill_manual(values = fill_colors)
  p
}

#' Waterfall line plot: cells retained at each filtering step
#'
#' @param counts_df  A tibble with columns `sample`, `step`, `n_cells`.
#' @param cfg        Config list from [sc_config()].
#' @return A ggplot object.
#' @export
sc_plot_filter_waterfall <- function(counts_df, cfg) {
  sample_colors <- cfg$.colors$samples

  ggplot2::ggplot(counts_df,
                  ggplot2::aes(x = step, y = n_cells,
                               group = sample, color = sample)) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::scale_color_manual(values = sample_colors) +
    ggplot2::labs(title = "Cells retained at each filtering step",
                  x = "Step", y = "Cell count") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x    = ggplot2::element_text(angle = 30, hjust = 1),
                   legend.position = "right")
}
