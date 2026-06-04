# =============================================================================
# R/preprocess.R
# Filtering, ambient RNA decontamination, doublet detection, merge.
# =============================================================================

#' Apply a permissive pre-filter to remove true empty droplets
#'
#' Retains cells with at least 100 genes and 200 UMIs. This low threshold
#' avoids biasing ambient RNA estimation by decontX.
#'
#' @param obj  A Seurat object.
#' @return Filtered Seurat object.
#' @export
sc_prefilter <- function(obj) {
  keep <- obj$nFeature_RNA >= 100L & obj$nCount_RNA >= 200L
  obj[, keep]
}

#' Run decontX ambient RNA decontamination on all samples
#'
#' When `raw_path` is specified in `params.yml`, empty droplets are loaded
#' and used as background, significantly improving estimation accuracy.
#'
#' If `decontX$use_decontX_counts: true` in `params.yml`, the decontaminated
#' count matrix is stored as a `"decontX"` assay in each object and used
#' as the default assay for downstream normalisation.
#'
#' @param obj_list  Named list of pre-filtered Seurat objects.
#' @param cfg       Config list from [sc_config()].
#' @param force     If `TRUE`, ignore cache.
#'
#' @return Named list of Seurat objects with `decontX_contamination` and
#'   `decontX_background` metadata columns added.
#' @export
sc_decontX <- function(obj_list, cfg, force = FALSE) {

  if (!isTRUE(cfg$decontX$run)) {
    cli::cli_alert_info("decontX disabled in params.yml — skipping")
    for (nm in names(obj_list)) {
      obj_list[[nm]]$decontX_contamination <- NA_real_
      obj_list[[nm]]$decontX_background    <- FALSE
    }
    return(obj_list)
  }

  sc_cache_run("02a_decontX_objects.qs2",
               cache_dir = sc_cache_dir(cfg),
               force     = force,
               expr = {
                 cli::cli_h1("Running decontX")
                 meta <- sc_sample_meta(cfg)

                 obj_list <- lapply(names(obj_list), function(sname) {
                   obj         <- obj_list[[sname]]
                   sample_info <- meta[meta$name == sname, ]
                   cli::cli_alert_info("{sname} ({ncol(obj)} cells)")
                   .run_decontX_one(obj, sample_info, cfg)
                 }) |> stats::setNames(names(obj_list))

                 n_bg <- sum(sapply(obj_list, function(o) {
                   isTRUE(unique(o$decontX_background))
                 }))
                 use_counts <- isTRUE(cfg$decontX$use_decontX_counts)
                 cli::cli_alert_success(
                   "decontX done — {n_bg}/{length(obj_list)} samples ran \\
                    with raw background | decontX counts used: {use_counts}"
                 )
                 obj_list
               })
}

#' Run scDblFinder doublet detection on all samples
#'
#' @param obj_list  Named list of Seurat objects.
#' @param cfg       Config list from [sc_config()].
#' @param force     If `TRUE`, ignore cache.
#'
#' @return Named list with `scDblFinder_score` and `scDblFinder_class` added.
#' @export
sc_doublets <- function(obj_list, cfg, force = FALSE) {

  sc_cache_run("02b_doublet_objects.qs2",
               cache_dir = sc_cache_dir(cfg),
               force     = force,
               expr = {
                 cli::cli_h1("Running scDblFinder")

                 obj_list <- lapply(names(obj_list), function(sname) {
                   obj <- obj_list[[sname]]
                   cli::cli_alert_info("{sname} ({ncol(obj)} cells)")

                   # Always use raw RNA counts for doublet detection
                   counts_mat <- Seurat::GetAssayData(
                     obj, assay = "RNA", layer = "counts"
                   )

                   sce <- scDblFinder::scDblFinder(
                     SingleCellExperiment::SingleCellExperiment(
                       list(counts = counts_mat)
                     ),
                     dbr.sd  = 0.015,
                     BPPARAM = BiocParallel::SerialParam(),
                     verbose = FALSE
                   )

                   obj$scDblFinder_score <- sce$scDblFinder.score
                   obj$scDblFinder_class <- sce$scDblFinder.class
                   obj
                 }) |> stats::setNames(names(obj_list))

                 n_dbl <- sapply(obj_list, function(o) {
                   sum(o$scDblFinder_class == "doublet")
                 })
                 cli::cli_alert_success(
                   "Doublets: {paste(names(n_dbl), n_dbl, sep = '=', collapse = ', ')}"
                 )
                 obj_list
               })
}

#' Apply final QC filters
#'
#' Applies per-sample thresholds for nFeature, nCount, %MT, %RB, complexity,
#' decontX contamination score, and optionally removes doublets.
#'
#' @param obj_list         Named list of Seurat objects (after doublet detection).
#' @param thresholds_list  Output of [sc_qc_all_thresholds()].
#' @param remove_doublets  If `TRUE`, remove cells classified as doublets.
#'
#' @return Named list of filtered Seurat objects.
#' @export
sc_filter <- function(obj_list, thresholds_list, remove_doublets = TRUE) {

  cli::cli_h1("Applying QC filters")

  lapply(names(obj_list), function(sname) {
    obj <- obj_list[[sname]]
    thr <- thresholds_list[[sname]]
    n0  <- ncol(obj)

    keep <- obj$nFeature_RNA >= thr$min_features &
      obj$nFeature_RNA <= thr$max_features &
      obj$nCount_RNA   >= thr$min_counts   &
      obj$nCount_RNA   <= thr$max_counts   &
      obj$log10_genes_per_umi >= thr$min_log10_genes_per_umi &
      (is.na(obj$decontX_contamination) |
         obj$decontX_contamination <= (thr$max_conta %||% 0.75))
    
    # Apply MT filter only if column exists and has signal
    if ("percent.mt" %in% colnames(obj@meta.data) &&
        !all(is.na(obj$percent.mt)) &&
        !all(obj$percent.mt == 0, na.rm = TRUE)) {
      keep <- keep & obj$percent.mt <= thr$max_mt_percent
    }
    
    # Apply RB filter only if column exists and has signal
    if ("percent.rb" %in% colnames(obj@meta.data) &&
        !all(is.na(obj$percent.rb)) &&
        !all(obj$percent.rb == 0, na.rm = TRUE)) {
      keep <- keep & obj$percent.rb <= thr$max_rb_percent
    }
    

    if (remove_doublets) keep <- keep & obj$scDblFinder_class == "singlet"

    n1 <- sum(keep)
    cli::cli_alert_info(
      "{sname}: {n0} \u2192 {n1} cells ({n0 - n1} removed, \\
       {round(100 * (n0 - n1) / n0, 1)}%)"
    )
    obj[, keep]
  }) |> stats::setNames(names(obj_list))
}

#' Merge all filtered objects into a single Seurat v5 object
#'
#' RNA layers remain separate after merge — [sc_integrate()] calls
#' `JoinLayers()` after integration.
#'
#' @param obj_list  Named list of filtered Seurat objects.
#' @param cfg       Config list from [sc_config()].
#' @param force     If `TRUE`, ignore cache.
#'
#' @return A merged Seurat v5 object.
#' @export
sc_merge <- function(obj_list, cfg, force = FALSE) {

  sc_cache_run("02c_merged.qs2",
               cache_dir = sc_cache_dir(cfg),
               force     = force,
               expr = {
                 cli::cli_h1("Merging {length(obj_list)} objects")

                 merged <- merge(
                   x            = obj_list[[1]],
                   y            = obj_list[-1],
                   add.cell.ids = names(obj_list)
                 )

                 cli::cli_alert_success(
                   "Merged: {ncol(merged)} cells | {nrow(merged)} genes | \\
                    {length(SeuratObject::Layers(merged, 'counts'))} RNA layers"
                 )
                 merged
               })
}

# ── Internal ──────────────────────────────────────────────────────────────────

#' @keywords internal
.run_decontX_one <- function(obj, sample_info, cfg) {

  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = Seurat::GetAssayData(obj, assay = "RNA", layer = "counts"))
  )

  # Remove all-zero genes — decontX produces NA on them
  nz <- Matrix::rowSums(SingleCellExperiment::counts(sce)) > 0
  if (any(!nz)) {
    cli::cli_alert_info("  Removing {sum(!nz)} all-zero genes before decontX")
    sce <- sce[nz, ]
  }

  bg_mat <- sc_load_background(sample_info,
                                SingleCellExperiment::counts(sce))

  result <- tryCatch({
    if (!is.null(bg_mat)) {
      common  <- intersect(rownames(sce), rownames(bg_mat))
      sce_aln <- sce[common, ]
      sce_bg  <- SingleCellExperiment::SingleCellExperiment(
        list(counts = bg_mat[common, , drop = FALSE])
      )
      sce_aln <- celda::decontX(sce_aln, background = sce_bg,
                                 delta = c(10, 10), verbose = FALSE)
      list(contamination = sce_aln$decontX_contamination,
           decontX_sce   = sce_aln,
           used_bg       = TRUE)
    } else {
      sce <- celda::decontX(sce, delta = c(10, 10), verbose = FALSE)
      list(contamination = sce$decontX_contamination,
           decontX_sce   = sce,
           used_bg       = FALSE)
    }
  }, error = function(e) {
    cli::cli_alert_warning(
      "decontX failed for {sample_info$name}: {conditionMessage(e)}"
    )
    list(contamination = rep(NA_real_, ncol(obj)),
         decontX_sce   = NULL,
         used_bg       = FALSE)
  })

  obj$decontX_contamination <- result$contamination
  obj$decontX_background    <- result$used_bg

  # Store decontaminated counts as a separate assay if requested
  if (isTRUE(cfg$decontX$use_decontX_counts) && !is.null(result$decontX_sce)) {
    decontX_counts <- celda::decontXcounts(result$decontX_sce)

    # Restore full gene set (decontX may have subset genes)
    full_counts <- Matrix::Matrix(0,
                                   nrow = nrow(obj),
                                   ncol = ncol(obj),
                                   dimnames = list(rownames(obj), colnames(obj)),
                                   sparse = TRUE)
    shared <- intersect(rownames(decontX_counts), rownames(obj))
    full_counts[shared, ] <- decontX_counts[shared, ]

    obj[["decontX"]] <- SeuratObject::CreateAssay5Object(
      counts = round(full_counts)
    )
    cli::cli_alert_info("  decontX counts stored as 'decontX' assay")
  }

  obj
}
