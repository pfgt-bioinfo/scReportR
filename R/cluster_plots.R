# =============================================================================
# R/cluster_plots.R
# Clustering and marker visualisation — all return ggplot objects.
# =============================================================================

# ── HVG plot ─────────────────────────────────────

#' Highly variable genes plot
#'
#' Handles both SCTransform v2 and LogNormalize workflows automatically.
#'
#' @param merged  Normalised Seurat object.
#' @param n_label Number of top HVGs to label.
#' @return A ggplot object.
#' @export
sc_plot_hvg <- function(merged, n_label = 20L) {
  
  default_assay <- Seurat::DefaultAssay(merged)
  
  method <- dplyr::case_when(
    "SCT"     %in% names(merged@assays) ~ "SCT",
    default_assay == "decontX"          ~ "decontX",
    TRUE                                ~ "LogNormalize"
  )
  
  if (method == "SCT") {
    hvg      <- Seurat::VariableFeatures(merged)
    attrs    <- merged[["SCT"]]@SCTModel.list[[1]]@feature.attributes
    hvg_df   <- data.frame(
      gene              = rownames(attrs),
      mean              = attrs$gmean,
      residual_variance = attrs$residual_variance,
      is_hvg            = rownames(attrs) %in% hvg
    ) |> dplyr::arrange(dplyr::desc(residual_variance))
    
    top_genes <- head(hvg_df$gene[hvg_df$is_hvg], n_label)
    
    ggplot2::ggplot(hvg_df,
                    ggplot2::aes(x     = mean,
                                 y     = residual_variance,
                                 color = is_hvg,
                                 label = gene)) +
      ggplot2::geom_point(size = 0.6, alpha = 0.6) +
      ggrepel::geom_text_repel(
        data         = dplyr::filter(hvg_df, gene %in% top_genes),
        size         = 3,
        max.overlaps = 20
      ) +
      ggplot2::scale_color_manual(
        values = c("FALSE" = "black", "TRUE" = "red"),
        labels = c("Other", "HVG"),
        name   = NULL
      ) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(title = "Highly variable genes (SCTransform v2)",
                    x     = "Geometric mean expression",
                    y     = "Residual variance") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "bottom")
    
  } else {
    # LogNormalize — use scCustomize
    meta  <- Seurat::HVFInfo(merged, assay = default_assay)
    hvg   <- Seurat::VariableFeatures(merged)  # non triés en v5
    
    # Trier explicitement par variance standardisée décroissante
    top_genes <- meta[hvg, ] |>
      dplyr::arrange(dplyr::desc(variance.standardized)) |>
      head(n_label) |>
      rownames()
    
    scCustomize::VariableFeaturePlot_scCustom(
      seurat_object = merged,
      num_features  = n_label,
      repel         = TRUE
    )
  }
}

# ── Generic dimensionality reduction plot ─────────────────────────────────────

#' Plot cells in a dimensionality reduction coloured by metadata or gene expression
#'
#' Generic function underlying [sc_plot_umap_cols()] and [sc_plot_pca_cols()].
#' Handles three cases automatically:
#' - **Categorical metadata** → [scCustomize::DimPlot_scCustom()]
#' - **Continuous metadata** (QC metrics, scores) → [scCustomize::FeaturePlot_scCustom()]
#' - **Gene names** → [scCustomize::FeaturePlot_scCustom()]
#'
#' @param merged     Seurat object.
#' @param cols       Character vector of metadata column names or gene names.
#' @param reduction  Dimensionality reduction to use (e.g. `"umap"`, `"pca"`).
#' @param cfg        Config list from [sc_config()]. Used for colour palettes.
#' @param pt_size    Point size.
#' @return Named list of ggplot objects, one per column/gene.
#' @export
sc_plot_dim_cols <- function(merged, cols, reduction = "umap",
                              cfg = NULL, pt_size = 0.3) {

  meta_cols <- cols[cols %in% colnames(merged@meta.data)]
  gene_cols <- cols[cols %in% rownames(merged)]
  missing   <- setdiff(cols, c(meta_cols, gene_cols))

  if (length(missing) > 0)
    cli::cli_alert_warning(
      "Not found in metadata or features: {paste(missing, collapse = ', ')}"
    )

  # Metadata plots
  meta_plots <- lapply(meta_cols, function(col) {
    is_continuous <- is.numeric(merged@meta.data[[col]])

    if (is_continuous) {
      scCustomize::FeaturePlot_scCustom(
        seurat_object = merged,
        features      = col,
        reduction     = reduction,
        pt.size       = pt_size,
        raster        = TRUE
      )
    } else {
      colors <- .resolve_palette(col, cfg)
      scCustomize::DimPlot_scCustom(
        seurat_object = merged,
        reduction     = reduction,
        group.by      = col,
        label         = TRUE,
        label.size    = 3,
        repel         = TRUE,
        raster        = TRUE,
        pt.size       = pt_size,
        colors_use    = colors,
        figure_plot   = TRUE
      )
    }
  }) |> stats::setNames(meta_cols)

  # Gene expression plots — always FeaturePlot
  gene_plots <- lapply(gene_cols, function(col) {
    scCustomize::FeaturePlot_scCustom(
      seurat_object = merged,
      features      = col,
      reduction     = reduction,
      pt.size       = pt_size,
      raster        = TRUE
    )
  }) |> stats::setNames(gene_cols)

  c(meta_plots, gene_plots)
}

#' Plot UMAP coloured by metadata or gene expression
#'
#' Wrapper around [sc_plot_dim_cols()] with `reduction = "umap"`.
#'
#' @inheritParams sc_plot_dim_cols
#' @export
sc_plot_umap_cols <- function(merged, cols, cfg = NULL, pt_size = 0.3) {
  sc_plot_dim_cols(merged, cols, reduction = "umap", cfg = cfg, pt_size = pt_size)
}

#' Plot PCA coloured by metadata or gene expression
#'
#' Wrapper around [sc_plot_dim_cols()] with `reduction = "pca"`.
#'
#' @inheritParams sc_plot_dim_cols
#' @export
sc_plot_pca_cols <- function(merged, cols, cfg = NULL, pt_size = 0.3) {
  sc_plot_dim_cols(merged, cols, reduction = "pca", cfg = cfg, pt_size = pt_size)
}

# ── Tabset helpers ─────────────────────────────────────────────────────────────

#' Render a named list of plots as a Quarto panel-tabset
#'
#' Must be called inside a chunk with `results='asis'`.
#'
#' @param plot_list  Named list of ggplot objects.
#' @param print_fn   Function used to render each plot (default: `print`).
#' @export
sc_tabset <- function(plot_list, print_fn = print) {
  cat("\n::: {.panel-tabset}\n")
  for (nm in names(plot_list)) {
    cat("\n##", nm, "\n\n")
    print_fn(plot_list[[nm]])
    cat("\n")
  }
  cat("\n:::\n")
}

#' UMAP tabset for all clustering resolutions
#'
#' @param merged       Clustered Seurat object.
#' @param resolutions  Numeric vector of resolutions to display.
#' @return Renders a tabset directly — call inside a chunk with `results='asis'`.
#' @export
sc_tabset_resolutions <- function(merged, resolutions) {
  cols <- paste0("clusters_res", resolutions)
  cols <- cols[cols %in% colnames(merged@meta.data)]

  plots <- lapply(cols, function(col) {
    n   <- length(unique(merged@meta.data[[col]]))
    res <- gsub("clusters_res", "", col)

    scCustomize::DimPlot_scCustom(
      seurat_object = merged,
      group.by      = col,
      label         = TRUE,
      label.size    = 3,
      repel         = TRUE,
      raster        = TRUE,
      pt.size       = 0.3,
      figure_plot   = TRUE
    ) +
      ggplot2::ggtitle(paste0("Res = ", res, " (", n, " clusters)")) +
      ggplot2::theme(legend.position = "none")
  }) |> stats::setNames(paste0("Res ", resolutions))

  sc_tabset(plots)
}

# ── Elbow plot ─────────────────────────────────────────────────────────────────

#' Elbow plot for PC selection
#'
#' @param merged  Seurat object after PCA.
#' @param cfg     Config list from [sc_config()].
#' @param n_show  Number of PCs to display.
#' @return A ggplot object.
#' @export
sc_plot_elbow <- function(merged, cfg, n_show = 50L) {
  pcs_use <- cfg$normalization$pcs_use %||% 30L

  Seurat::ElbowPlot(merged,
                    ndims = min(n_show, length(merged@reductions$pca))) +
    ggplot2::labs(title = "Elbow plot — PC selection") +
    ggplot2::geom_vline(xintercept = pcs_use,
                        linetype = "dashed", color = "red") +
    ggplot2::theme_bw()
}

# ── Cluster composition ────────────────────────────────────────────────────────

#' Stacked bar chart of cluster composition by metadata variables
#'
#' @param merged      Clustered Seurat object.
#' @param cols        Character vector of metadata columns to plot.
#' @param cfg         Config list from [sc_config()].
#' @param reverse     If `TRUE`, metadata on x axis and clusters as fill.
#' @param legend_pos  Legend position.
#' @return Renders a tabset directly — call inside a chunk with `results='asis'`.
#' @export
sc_plot_cluster_composition <- function(
    merged, cols, cfg,
    reverse    = FALSE,
    legend_pos = c("bottom", "top", "left", "right", "none")) {

  cluster_col    <- paste0("clusters_res",
                            cfg$clustering$default_resolution %||% 0.6)
  legend_pos     <- match.arg(legend_pos)
  cluster_levels <- sort(unique(merged@meta.data[[cluster_col]]))
  cluster_colors <- stats::setNames(scales::hue_pal()(length(cluster_levels)),
                                    cluster_levels)

  # Filter to columns actually present in the object
  cols <- cols[cols %in% colnames(merged@meta.data)]
  if (length(cols) == 0L) {
    cli::cli_alert_warning("None of the requested columns found in metadata")
    return(invisible(NULL))
  }

  plots <- lapply(cols, function(col) {
    meta_colors <- .resolve_palette(col, cfg)

    if (!reverse) {
      df <- merged@meta.data |>
        dplyr::count(.data[[cluster_col]], .data[[col]]) |>
        dplyr::group_by(.data[[cluster_col]]) |>
        dplyr::mutate(pct = 100 * n / sum(n))

      p <- ggplot2::ggplot(df,
                           ggplot2::aes(x    = .data[[cluster_col]],
                                        y    = pct,
                                        fill = .data[[col]])) +
        ggplot2::labs(title = paste("Cluster composition by", col),
                      x = "Cluster", y = "% cells", fill = col)
      if (!is.null(meta_colors))
        p <- p + ggplot2::scale_fill_manual(values = meta_colors)

    } else {
      df <- merged@meta.data |>
        dplyr::count(.data[[col]], .data[[cluster_col]]) |>
        dplyr::group_by(.data[[col]]) |>
        dplyr::mutate(pct = 100 * n / sum(n))

      p <- ggplot2::ggplot(df,
                           ggplot2::aes(x    = .data[[col]],
                                        y    = pct,
                                        fill = .data[[cluster_col]])) +
        ggplot2::scale_fill_manual(values = cluster_colors) +
        ggplot2::labs(title = paste("Composition of", col, "by cluster"),
                      x = col, y = "% cells", fill = "Cluster")
    }

    p +
      ggplot2::geom_col(alpha = 0.85) +
      ggplot2::scale_y_continuous(labels = function(x) paste0(x, "%")) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x     = ggplot2::element_text(angle = 45,
                                                              hjust = 1),
                     legend.position = legend_pos)
  }) |> stats::setNames(cols)

  sc_tabset(plots)
}

# ── Marker plots ───────────────────────────────────────────────────────────────

#' DotPlot of markers across clusters
#'
#' Can display either known canonical markers (from `params.yml`) or
#' differential markers computed by [sc_find_markers()].
#'
#' @param merged       Seurat object.
#' @param cfg          Config list from [sc_config()].
#' @param markers_list Named list of character vectors (known markers).
#'   If provided, `markers_df` is ignored.
#' @param markers_df   Output of [sc_find_markers()]. Used when
#'   `markers_list` is `NULL`.
#' @param n_top        Number of top markers per cluster to display
#'   when using `markers_df`.
#' @return A plot object.
#' @export
sc_plot_markers_dot <- function(merged, cfg, markers_list = NULL,
                                markers_df = NULL, n_top = 5L) {
  
  col <- paste0("clusters_res", cfg$clustering$default_resolution %||% 0.6)
  
  if (!is.null(markers_list)) {
    # Known canonical markers from params.yml
    genes <- unique(unlist(markers_list))
    
  } else if (!is.null(markers_df) && nrow(markers_df) > 0L &&
             "cluster" %in% colnames(markers_df)) {
    # Computed differential markers
    genes <- markers_df |>
      dplyr::group_by(cluster) |>
      dplyr::slice_max(order_by = avg_log2FC, n = n_top, with_ties = FALSE) |>
      dplyr::pull(gene) |>
      unique()
    
  } else {
    cli::cli_alert_warning(
      "Provide either markers_list or markers_df to sc_plot_markers_dot()"
    )
    return(NULL)
  }
  
  genes <- genes[genes %in% rownames(merged)]
  
  scCustomize::Clustered_DotPlot(
    seurat_object = merged,
    features      = genes,
    group.by      = col,
    plot_km_elbow = FALSE
  )
}

#' Heatmap of top differential markers per cluster
#'
#' Works with both SCT and RNA assays. When SCT is used, only genes present
#' in `scale.data` are displayed; remaining genes are silently dropped.
#'
#' @param merged      Seurat object.
#' @param markers_df  Output of [sc_find_markers()].
#' @param n_top       Number of top markers per cluster.
#' @return A ggplot object, or `NULL` if no markers available.
#' @export
sc_plot_markers_heatmap <- function(merged, markers_df, n_top = 5L) {
  
  if (is.null(markers_df) || nrow(markers_df) == 0L ||
      !"cluster" %in% colnames(markers_df)) {
    cli::cli_alert_warning("No markers to plot")
    return(NULL)
  }
  
  top <- markers_df |>
    dplyr::group_by(cluster) |>
    dplyr::slice_max(order_by = avg_log2FC, n = n_top, with_ties = FALSE) |>
    dplyr::pull(gene) |>
    unique()
  
  # Use the assay that has scale.data available
  # Priority: RNA (always scaled) > decontX > SCT
  assay_use <- dplyr::case_when(
    "scale.data" %in% SeuratObject::Layers(merged[["RNA"]])    ~ "RNA",
    "decontX" %in% names(merged@assays) &&
      "scale.data" %in% SeuratObject::Layers(merged[["decontX"]]) ~ "decontX",
    TRUE ~ "RNA"
  )
  
  cli::cli_alert_info("Heatmap using assay: {assay_use}")
  prev_assay <- Seurat::DefaultAssay(merged)
  Seurat::DefaultAssay(merged) <- assay_use
  
  # Filter to genes present in scale.data
  scaled_genes <- rownames(merged[[assay_use]]@scale.data)
  missing      <- top[!top %in% scaled_genes]
  if (length(missing) > 0L)
    cli::cli_alert_info("{length(missing)} markers not in scale.data — omitted")
  top <- top[top %in% scaled_genes]
  
  if (length(top) == 0L) {
    cli::cli_alert_warning("No markers remain after scale.data filtering")
    Seurat::DefaultAssay(merged) <- prev_assay
    return(NULL)
  }
  
  cells_use <- scCustomize::Random_Cells_Downsample(
    seurat_object = merged, num_cells = 200, allow_lower = TRUE
  )
  
  p <- Seurat::DoHeatmap(merged, features = top, cells = cells_use,
                         slot = "scale.data", size = 3L, angle = 90) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6)) +
    ggplot2::labs(title = paste0("Top ", n_top, " markers per cluster"))
  
  Seurat::DefaultAssay(merged) <- prev_assay
  p
}

#' Clustree plot for resolution stability
#'
#' @param merged       Clustered Seurat object.
#' @param resolutions  Numeric vector of resolutions tested.
#' @return A ggplot object, or `NULL` if clustree is not installed.
#' @export
sc_plot_clustree <- function(merged, resolutions) {
  if (!requireNamespace("clustree", quietly = TRUE)) {
    cli::cli_alert_warning("Package 'clustree' not installed")
    return(NULL)
  }
  clustree::clustree(merged@meta.data,
                     prefix      = "clusters_res",
                     edge_colour = "grey60",
                     edge_arrow  = FALSE) +
    ggplot2::labs(title = "Clustree — cluster stability across resolutions")
}

#' UMAP tabset of module scores for each cluster
#'
#' @param merged  Seurat object with module scores from [sc_module_scores()].
#' @param cfg     Config list from [sc_config()].
#' @return Renders a tabset directly — call inside a chunk with `results='asis'`.
#' @export
sc_plot_module_scores <- function(merged, cfg) {

  score_cols <- grep("^score_cluster_", colnames(merged@meta.data), value = TRUE)

  if (length(score_cols) == 0L) {
    cli::cli_alert_warning(
      "No module score columns found — run sc_module_scores() first"
    )
    return(invisible(NULL))
  }

  score_cols <- score_cols[order(
    as.numeric(gsub("score_cluster_", "", score_cols))
  )]

  sc_plot_dim_cols(merged, cols = score_cols) |>
    stats::setNames(gsub("score_cluster_", "Cluster ", score_cols)) |>
    sc_tabset()
}
