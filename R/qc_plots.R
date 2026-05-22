# =============================================================================
# R/qc_plots.R
# QC visualisation functions — all return ggplot objects.
# =============================================================================

#' Generic bar plot of cell counts with configurable fill
#'
#' Accepts a Seurat object, a named list of Seurat objects, or a plain
#' data frame. The `fill` column drives colours; pass `cfg` to use the
#' project palettes automatically.
#'
#' @param obj_or_meta  A Seurat object, named list of Seurat objects, or data frame.
#' @param x            Column name for the x axis.
#' @param fill         Column name for the fill colour.
#' @param facet        Optional column name to facet by.
#' @param cfg          Config list from [sc_config()]. Used for colour palettes.
#' @param title        Plot title.
#' @param position     Bar position: `"dodge"` or `"stack"`. Auto-detected
#'   from `facet` when `NULL`.
#' @return A ggplot object.
#' @export
sc_plot_cell_counts <- function(obj_or_meta, x, fill, facet = NULL,
                                cfg = NULL, title = "Cell counts",
                                position = NULL) {

  meta <- if (inherits(obj_or_meta, "Seurat")) {
    as.data.frame(obj_or_meta@meta.data)
  } else if (is.list(obj_or_meta) && !inherits(obj_or_meta, "Seurat")) {
    dplyr::bind_rows(lapply(obj_or_meta, function(o) as.data.frame(o@meta.data)))
  } else {
    obj_or_meta
  }

  counts <- meta |>
    dplyr::count(dplyr::across(dplyr::all_of(unique(c(x, fill, facet)))))

  pos        <- position %||% if (is.null(facet)) "dodge" else "stack"
  pos_geom   <- if (pos == "dodge") ggplot2::position_dodge(width = 0.75)
                else ggplot2::position_stack()
  pos_text   <- if (pos == "dodge") ggplot2::position_dodge(width = 0.75)
                else ggplot2::position_stack(vjust = 1.02)
  text_vjust <- if (pos == "dodge") -0.4 else 0

  colors <- .resolve_palette(fill, cfg)

  p <- ggplot2::ggplot(counts,
                       ggplot2::aes(x = .data[[x]], y = n,
                                    fill = .data[[fill]])) +
    ggplot2::geom_col(alpha = 0.85, width = 0.7, position = pos_geom) +
    ggplot2::geom_text(ggplot2::aes(label = scales::comma(n)),
                       position = pos_text, vjust = text_vjust, size = 3) +
    ggplot2::scale_y_continuous(labels = scales::comma,
                                expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(title = title, x = NULL, y = "Cell count", fill = fill) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x     = ggplot2::element_text(angle = 45,
                                                            hjust = 1, size = 8),
                   legend.position = "bottom")

  if (!is.null(colors)) p <- p + ggplot2::scale_fill_manual(values = colors)

  if (!is.null(facet)) {
    p <- p +
      ggplot2::facet_grid(stats::as.formula(paste("~", facet)),
                          scales = "free_x", space = "free_x") +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill = "grey90"))
  }
  p
}

#' Strip plot of median UMI and gene counts per sample
#'
#' @param obj_list  Named list of Seurat objects.
#' @param cfg       Config list from [sc_config()].
#' @return A ggplot object.
#' @export
sc_plot_qc_distri <- function(obj_list, cfg) {
  sample_colors <- cfg$.colors$samples

  sc_qc_summary(obj_list) |>
    tidyr::pivot_longer(cols      = c(median_umi, median_genes),
                        names_to  = "metric",
                        values_to = "value") |>
    ggplot2::ggplot(ggplot2::aes(x     = metric,
                                 y     = value,
                                 color = sample,
                                 label = sample)) +
    ggplot2::geom_jitter(size = 3, width = 0.1, alpha = 0.8) +
    ggplot2::scale_y_continuous(labels = scales::comma) +
    ggplot2::scale_color_manual(values = sample_colors) +
    ggplot2::facet_wrap(~metric, scales = "free_y") +
    ggplot2::labs(x = NULL, y = NULL,
                  title = "QC metrics distribution across samples") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x     = ggplot2::element_blank(),
                   axis.ticks.x    = ggplot2::element_blank(),
                   legend.position = "right")
}

#' Combined QC violin plot — all samples, 4 panels
#'
#' Threshold lines are drawn as per-sample segments so that per-sample
#' overrides are immediately visible.
#'
#' @param obj_list         Named list of Seurat objects.
#' @param thresholds_list  Output of [sc_qc_all_thresholds()].
#' @param cfg              Config list from [sc_config()].
#' @param fill_by          Metadata column used for fill colour. Defaults to
#'   the best available categorical column (library_type > tissue > population).
#' @return A patchwork ggplot.
#' @export
sc_plot_qc_violins <- function(obj_list, thresholds_list, cfg,
                                fill_by = NULL) {

  sample_levels <- names(obj_list)
  meta_all <- dplyr::bind_rows(
    lapply(obj_list, function(o) as.data.frame(o@meta.data))
  ) |>
    dplyr::mutate(sample_id = factor(sample_id, levels = sample_levels))

  # Auto-detect fill column if not provided
  fill_col <- fill_by %||% .best_meta_col(meta_all)
  if (is.null(fill_col)) {
    meta_all$.fill_dummy <- "all"
    fill_col <- ".fill_dummy"
  }

  fill_colors <- .resolve_palette(fill_col, cfg)

  # Per-sample threshold segments (x position = integer factor level)
  thr_df <- dplyr::bind_rows(lapply(seq_along(sample_levels), function(i) {
    thr <- thresholds_list[[sample_levels[i]]]
    tibble::tibble(
      sample_id    = sample_levels[i],
      x_num        = i,
      min_features = thr$min_features,
      max_features = thr$max_features,
      min_counts   = thr$min_counts,
      max_counts   = thr$max_counts,
      max_mt       = thr$max_mt_percent,
      max_rb       = thr$max_rb_percent
    )
  }))

  .make_panel <- function(yvar, ylab, thr_cols) {
    p <- ggplot2::ggplot(
      meta_all,
      ggplot2::aes(x = sample_id, y = .data[[yvar]],
                   fill = .data[[fill_col]])
    ) +
      ggplot2::geom_violin(scale = "width", alpha = 0.75, trim = TRUE,
                           linewidth = 0.3) +
      ggplot2::labs(x = NULL, y = ylab, title = ylab) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
        legend.position = "none",
        plot.title      = ggplot2::element_text(size = 9, face = "bold")
      )

    if (!is.null(fill_colors))
      p <- p + ggplot2::scale_fill_manual(values = fill_colors)

    for (col in thr_cols) {
      seg <- thr_df |> dplyr::select(x_num, y = dplyr::all_of(col))
      p <- p + ggplot2::geom_segment(
        data        = seg,
        ggplot2::aes(x = x_num - 0.38, xend = x_num + 0.38, y = y, yend = y),
        color       = "red", linetype = "dashed", linewidth = 0.55,
        inherit.aes = FALSE
      )
    }
    p
  }

  p1 <- .make_panel("nFeature_RNA", "Genes / cell",    c("min_features", "max_features"))
  p2 <- .make_panel("nCount_RNA",   "UMI / cell",      c("min_counts",   "max_counts"))
  p3 <- .make_panel("percent.mt",   "% Mitochondrial", "max_mt")
  p4 <- .make_panel("percent.rb",   "% Ribosomal",     "max_rb")

  (p1 + p2 + p3 + p4) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
}

#' Scatter plots: UMI vs genes, UMI vs MT%, and cell complexity
#'
#' Uses [scCustomize::QC_Plot_UMIvsGene()] for the main panel and adds
#' two complementary panels.
#'
#' @param obj_list         Named list of Seurat objects.
#' @param thresholds_list  Output of [sc_qc_all_thresholds()].
#' @return Named list of patchwork ggplots, one per sample.
#' @export
sc_plot_qc_scatter <- function(obj_list, thresholds_list) {

  lapply(names(obj_list), function(sname) {
    obj <- obj_list[[sname]]
    thr <- thresholds_list[[sname]]

    p1 <- scCustomize::QC_Plot_UMIvsGene(
      seurat_object      = obj,
      meta_gradient_name = "percent.mt",
      low_cutoff_gene    = thr$min_features,
      high_cutoff_gene   = thr$max_features,
      low_cutoff_UMI     = thr$min_counts,
      high_cutoff_UMI    = thr$max_counts
    ) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 11))

    p2 <- ggplot2::ggplot(as.data.frame(obj@meta.data),
                          ggplot2::aes(nCount_RNA, percent.mt)) +
      ggplot2::geom_point(size = 0.3, alpha = 0.4, color = "#2c7bb6") +
      ggplot2::geom_smooth(method = "loess", se = TRUE, color = "red",
                           linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = thr$max_mt_percent,
                          linetype = "dashed", color = "red") +
      ggplot2::scale_x_log10(labels = scales::comma) +
      ggplot2::labs(x = "UMI (log10)", y = "% Mitochondrial") +
      ggplot2::theme_bw()

    p3 <- ggplot2::ggplot(as.data.frame(obj@meta.data),
                          ggplot2::aes(nCount_RNA, log10_genes_per_umi)) +
      ggplot2::geom_point(size = 0.3, alpha = 0.4, color = "#5ab4ac") +
      ggplot2::geom_hline(yintercept = thr$min_log10_genes_per_umi,
                          linetype = "dashed", color = "red") +
      ggplot2::scale_x_log10(labels = scales::comma) +
      ggplot2::labs(x = "UMI (log10)",
                    y = "Complexity\n(log10 genes / log10 UMI)") +
      ggplot2::theme_bw()

    p1 + p2 + p3 +
      patchwork::plot_annotation(
        title = sname,
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold", size = 11)
        )
      )
  }) |> stats::setNames(names(obj_list))
}

#' Cell complexity density plot
#'
#' @param obj_list  Named list of Seurat objects.
#' @param cfg       Config list from [sc_config()].
#' @param fill_by   Metadata column used for fill. Defaults to best available.
#' @param facet_by  Metadata column used for faceting. Defaults to `"sample_id"`.
#' @return A ggplot object.
#' @export
sc_plot_complexity <- function(obj_list, cfg, fill_by = NULL,
                                facet_by = "sample_id") {

  meta_all <- dplyr::bind_rows(
    lapply(obj_list, function(o) as.data.frame(o@meta.data))
  )

  fill_col   <- fill_by %||% .best_meta_col(meta_all)
  fill_colors <- .resolve_palette(fill_col, cfg)

  # Validate facet column
  if (!facet_by %in% colnames(meta_all)) {
    cli::cli_alert_warning("Column '{facet_by}' not found — using sample_id")
    facet_by <- "sample_id"
  }

  p <- ggplot2::ggplot(meta_all,
                       ggplot2::aes(x = log10_genes_per_umi,
                                    fill = .data[[fill_col]])) +
    ggplot2::geom_density(alpha = 0.5) +
    ggplot2::geom_vline(xintercept = 0.8, linetype = "dashed", color = "red") +
    ggplot2::facet_wrap(stats::as.formula(paste("~", facet_by)),
                        scales = "free_y") +
    ggplot2::labs(title    = "Cell complexity",
                  subtitle = "Dashed = threshold 0.8",
                  x        = "log10(nGenes) / log10(nUMI)",
                  y        = "Density") +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "bottom")

  if (!is.null(fill_colors))
    p <- p + ggplot2::scale_fill_manual(values = fill_colors)
  p
}

#' Threshold overview plot (ranges per sample, overrides highlighted)
#'
#' @param thresholds_list  Output of [sc_qc_all_thresholds()].
#' @param cfg              Config list from [sc_config()].
#' @return A patchwork ggplot.
#' @export
sc_plot_thresholds <- function(thresholds_list, cfg) {
  thr_df <- sc_qc_thresholds_df(thresholds_list, cfg) |>
    dplyr::mutate(sample = factor(sample, levels = names(thresholds_list)))

  colors <- c("FALSE" = "grey50", "TRUE" = "#e74c3c")

  .bar <- function(ymin, ymax, title) {
    ggplot2::ggplot(thr_df,
                    ggplot2::aes(x = sample, color = has_override)) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = .data[[ymin]], ymax = .data[[ymax]]),
        width = 0.3, linewidth = 1
      ) +
      ggplot2::scale_color_manual(values = colors,
                                  labels = c("Default", "Override"),
                                  name   = NULL) +
      ggplot2::labs(title = title, x = NULL) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x     = ggplot2::element_text(angle = 45,
                                                              hjust = 1, size = 7),
                     legend.position = "bottom")
  }

  .dot <- function(y, title) {
    ggplot2::ggplot(thr_df,
                    ggplot2::aes(x = sample, y = .data[[y]],
                                 color = has_override)) +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_color_manual(values = colors, name = NULL) +
      ggplot2::labs(title = title, x = NULL) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x     = ggplot2::element_text(angle = 45,
                                                              hjust = 1, size = 7),
                     legend.position = "bottom")
  }

  (.bar("min_features", "max_features", "nFeature range") +
   .bar("min_counts",   "max_counts",   "nCount range")   +
   .dot("max_mt_percent", "Max MT %")) +
    patchwork::plot_annotation(
      title    = "Configured QC thresholds",
      subtitle = "Red = per-sample override in params.yml"
    )
}
