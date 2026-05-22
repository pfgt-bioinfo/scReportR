# =============================================================================
# R/annotate.R
# Automatic cell type annotation: SingleR + Pan-Human Azimuth + manual helper.
# =============================================================================

#' Run SingleR automatic annotation
#'
#' Works with both SCT and LogNormalize workflows. When SCT is the default
#' assay, log-normalised data is fetched from the RNA assay (layer `"data"`).
#'
#' @param merged  Seurat object (JoinLayers already applied).
#' @param cfg     Config list from [sc_config()].
#' @param force   If `TRUE`, ignore cache.
#'
#' @return Seurat object with SingleR labels added to metadata.
#' @export
sc_singler <- function(merged, cfg, force = FALSE) {

  if (!isTRUE(cfg$annotation$singler$run)) {
    cli::cli_alert_info("SingleR disabled in params.yml — skipping")
    return(merged)
  }

  sc_cache_run("05a_singler.qs2",
               cache_dir = sc_cache_dir(cfg),
               force     = force,
               expr = {
                 s <- cfg$annotation$singler
                 cli::cli_h1("SingleR")

                 # SingleR requires log-normalised counts.
                 # RNA assay "data" layer holds log-normalised data for both
                 # SCTransform (via NormalizeData called in sc_panhuman_azimuth
                 # or available from the original LogNormalize run) and
                 # LogNormalize workflows.
                 # If RNA "data" is empty, fall back to computing it on the fly.
                 rna_data <- tryCatch(
                   Seurat::GetAssayData(merged, assay = "RNA", layer = "data"),
                   error = function(e) NULL
                 )

                 if (is.null(rna_data) ||
                     all(rna_data@x == 0) ||
                     nrow(rna_data) == 0L) {
                   cli::cli_alert_info(
                     "RNA 'data' layer empty — running NormalizeData for SingleR"
                   )
                   prev_assay <- Seurat::DefaultAssay(merged)
                   Seurat::DefaultAssay(merged) <- "RNA"
                   merged     <- Seurat::NormalizeData(merged, verbose = FALSE)
                   Seurat::DefaultAssay(merged) <- prev_assay
                   rna_data   <- Seurat::GetAssayData(merged, assay = "RNA",
                                                       layer = "data")
                 }

                 for (ref_name in s$references) {
                   cli::cli_alert_info("Reference: {ref_name}")

                   ref <- do.call(
                     utils::getFromNamespace(ref_name, "celldex"),
                     list()
                   )

                   pred <- SingleR::SingleR(
                     test      = rna_data,
                     ref       = ref,
                     labels    = ref[[s$label_type %||% "label.main"]],
                     de.method = s$de_method %||% "wilcox",
                     BPPARAM   = BiocParallel::MulticoreParam(
                       min(4L, parallel::detectCores() - 1L)
                     )
                   )

                   prefix <- paste0("singler_",
                                    gsub("Data$", "", ref_name), "_")
                   merged@meta.data[[paste0(prefix, "label")]]  <- pred$labels
                   merged@meta.data[[paste0(prefix, "pruned")]] <- pred$pruned.labels
                   merged@meta.data[[paste0(prefix, "score")]]  <- apply(
                     pred$scores, 1L, max
                   )
                 }
                 cli::cli_alert_success(
                   "SingleR done — {length(s$references)} reference(s)"
                 )
                 merged
               })
}

#' Run Pan-Human Azimuth annotation via cloud API
#'
#' Requires the `AzimuthAPI` package and internet access. Data is processed on
#' Satijalab cloud servers. Always uses log-normalised RNA counts internally.
#'
#' @param merged  Seurat object.
#' @param cfg     Config list from [sc_config()].
#' @param force   If `TRUE`, ignore cache.
#'
#' @return Seurat object with Pan-Human Azimuth annotations in metadata.
#' @export
sc_panhuman_azimuth <- function(merged, cfg, force = FALSE) {

  if (!isTRUE(cfg$annotation$panhuman_azimuth$run)) {
    cli::cli_alert_info("Pan-Human Azimuth disabled in params.yml — skipping")
    return(merged)
  }

  if (!requireNamespace("AzimuthAPI", quietly = TRUE)) {
    cli::cli_abort(
      "Package 'AzimuthAPI' not installed.
       Install with: devtools::install_github('satijalab/AzimuthAPI')"
    )
  }

  sc_cache_run("05b_panhuman_azimuth.qs2",
               cache_dir = sc_cache_dir(cfg),
               force     = force,
               expr = {
                 cli::cli_h1("Pan-Human Azimuth (cloud API)")
                 cli::cli_alert_warning(
                   "Data will be sent to Satijalab cloud servers for annotation."
                 )

                 prev_assay <- Seurat::DefaultAssay(merged)

                 # Ensure log-normalised RNA data is available
                 if (prev_assay == "SCT") {
                   cli::cli_alert_info(
                     "SCT detected — running NormalizeData on RNA for Azimuth"
                   )
                   Seurat::DefaultAssay(merged) <- "RNA"
                   merged <- Seurat::NormalizeData(merged, verbose = FALSE)
                 }

                 merged <- AzimuthAPI::CloudAzimuth(merged)

                 # Restore default assay
                 Seurat::DefaultAssay(merged) <- prev_assay

                 az_cols <- grep("azimuth|full_hierarchical|softmax",
                                 colnames(merged@meta.data), value = TRUE)
                 cli::cli_alert_success(
                   "Pan-Human Azimuth done — {length(az_cols)} column(s) added"
                 )
                 merged
               })
}

#' Apply manual cluster-to-celltype annotation
#'
#' @param merged      Seurat object.
#' @param annotation  Named character vector: `c("1" = "T cells", ...)`.
#' @param cfg         Config list from [sc_config()].
#' @param new_col     Name of the new metadata column.
#'
#' @return Seurat object with manual annotation column added.
#' @export
sc_annotate_manual <- function(merged, annotation, cfg,
                               new_col = "celltype_manual") {
  cluster_col <- paste0("clusters_res",
                        cfg$clustering$default_resolution %||% 0.6)

  merged@meta.data[[new_col]] <- annotation[
    as.character(merged@meta.data[[cluster_col]])
  ]

  n <- sum(!is.na(merged@meta.data[[new_col]]))
  n_total <- length(unique(merged@meta.data[[cluster_col]]))
  cli::cli_alert_success("{n} / {n_total} clusters manually annotated")
  merged
}

#' Export final metadata to CSV
#'
#' Exports all QC, clustering, and annotation columns. Optional metadata
#' columns (tissue, population) are included only when present.
#'
#' @param merged  Annotated Seurat object.
#' @param cfg     Config list from [sc_config()].
#' @return Path to the exported CSV file (invisibly).
#' @export
sc_export_metadata <- function(merged, cfg) {
  path <- file.path(sc_output_dir(cfg), "metadata_final.csv")

  meta <- merged@meta.data |>
    tibble::rownames_to_column("cell_barcode") |>
    dplyr::select(
      cell_barcode,
      sample_id,
      dplyr::any_of(c("tissue", "library_type", "population", "species")),
      dplyr::starts_with("clusters_"),
      dplyr::contains("singler"),
      dplyr::contains("azimuth"),
      dplyr::any_of(c("celltype_manual",
                       "scDblFinder_class",
                       "scDblFinder_score",
                       "decontX_contamination",
                       "percent.mt",
                       "percent.rb",
                       "log10_genes_per_umi"))
    )

  utils::write.csv(meta, path, row.names = FALSE)
  cli::cli_alert_success(
    "Metadata exported: {nrow(meta)} cells \u2192 {.path {path}}"
  )
  invisible(path)
}
