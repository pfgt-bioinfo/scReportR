# =============================================================================
# R/normalize.R
# Normalisation (SCTransform v2 or LogNormalize) + PCA + integration.
# =============================================================================

#' Normalise a merged Seurat v5 object
#'
#' Applies SCTransform v2 (recommended) or LogNormalize + ScaleData,
#' then runs PCA. Results are cached.
#'
#' When `decontX$use_decontX_counts: true` in `params.yml` and the `decontX`
#' assay is present, it is used as the source of counts for normalisation
#' instead of the raw `RNA` assay.
#'
#' @param merged  Merged Seurat object from [sc_merge()].
#' @param cfg     Config list from [sc_config()].
#' @param force   If `TRUE`, ignore cache.
#'
#' @return Normalised Seurat object with PCA reduction.
#' @export
sc_normalize <- function(merged, cfg, force = FALSE) {
  
  sc_cache_run("03a_normalized.qs2",
               cache_dir = sc_cache_dir(cfg),
               force     = force,
               expr = {
                 p      <- cfg$normalization
                 method <- p$method %||% "SCTransform"
                 n_pcs  <- p$n_pcs  %||% 50L
                 
                 # Use decontX assay if requested and available
                 use_decontX  <- isTRUE(cfg$decontX$use_decontX_counts) &&
                   "decontX" %in% names(merged@assays)
                 source_assay <- if (use_decontX) "decontX" else "RNA"
                 
                 if (use_decontX)
                   cli::cli_alert_info("Using decontX counts for normalisation")
                 
                 # ── Step 1 : LogNormalize ──────────────────────────────────
                 # source_assay = "decontX" if use_decontX_counts, else "RNA"
                 cli::cli_h1("LogNormalize {source_assay} assay")
                 
                 merged <- Seurat::NormalizeData(merged, assay = source_assay,
                                                 verbose = FALSE)
                 merged <- Seurat::FindVariableFeatures(
                   merged,
                   assay     = source_assay,
                   nfeatures = p$n_variable_features %||% 3000L,
                   verbose   = FALSE
                 )
                 all_genes <- rownames(merged)
                 merged <- Seurat::ScaleData(
                   merged,
                   assay           = source_assay,
                   features        = all_genes,
                   vars.to.regress = unlist(p$vars_to_regress),
                   verbose         = FALSE
                 )
                 cli::cli_alert_success("{source_assay} normalised and scaled")
                 
                 # ── Step 2 : SCTransform (optional) ───────────────────────
                 if (method == "SCTransform") {
                   cli::cli_h1("SCTransform (assay = {source_assay})")
                   merged <- Seurat::SCTransform(
                     merged,
                     assay               = source_assay,
                     vst.flavor          = p$sct_flavor %||% "v2",
                     vars.to.regress     = unlist(p$vars_to_regress),
                     variable.features.n = p$n_variable_features %||% 3000L,
                     verbose             = FALSE
                   )
                   assay_pca <- "SCT"
                 } else {
                   assay_pca <- source_assay   # "decontX" ou "RNA"
                 }
                 
                 # ── Step 3 : PCA ───────────────────────────────────────────
                 cli::cli_alert_info("PCA: {n_pcs} PCs on {assay_pca}")
                 merged <- Seurat::RunPCA(merged, assay = assay_pca,
                                          npcs = n_pcs, verbose = FALSE)
                 
                 
                 # Set default assay for downstream analyses
                 default_assay <- dplyr::case_when(
                   "SCT"     %in% names(merged@assays) ~ "SCT",
                   "decontX" %in% names(merged@assays) &&
                     isTRUE(cfg$decontX$use_decontX_counts) ~ "decontX",
                   TRUE ~ "RNA"
                 )
                 Seurat::DefaultAssay(merged) <- default_assay
                 cli::cli_alert_info("DefaultAssay set to: {default_assay}")
                 
                 cli::cli_alert_success("Normalisation + PCA done")
                 merged
               })
}

#' Integrate samples and compute UMAP
#'
#' Supports Harmony (default), RPCA, CCA, or no integration.
#' Joins RNA layers after integration (required for FindMarkers).
#'
#' @param merged  Normalised Seurat object from [sc_normalize()].
#' @param cfg     Config list from [sc_config()].
#' @param force   If `TRUE`, ignore cache.
#'
#' @return Integrated Seurat object with UMAP reduction.
#' @export
sc_integrate <- function(merged, cfg, force = FALSE) {

  sc_cache_run("03b_integrated.qs2",
               cache_dir = sc_cache_dir(cfg),
               force     = force,
               expr = {
                 int    <- cfg$integration
                 method <- int$method %||% "Harmony"
                 pcs    <- cfg$normalization$pcs_use %||% 30L
                 cli::cli_h1("Integration: {method}")

                 if (method == "none") {
                   merged <- Seurat::RunUMAP(
                     merged, reduction = "pca", dims = 1:pcs, verbose = FALSE
                   )
                 } else {
                   method_fn <- switch(method,
                     Harmony         = Seurat::HarmonyIntegration,
                     RPCAIntegration = Seurat::RPCAIntegration,
                     CCAIntegration  = Seurat::CCAIntegration,
                     cli::cli_abort("Unknown integration method: {method}")
                   )

                   extra <- if (method == "Harmony") {
                     list(theta            = int$harmony_theta %||% 2,
                          max.iter.harmony = int$harmony_max_iter %||% 10L)
                   } else list()

                   merged <- do.call(
                     Seurat::IntegrateLayers,
                     c(list(object         = merged,
                            method         = method_fn,
                            orig.reduction = "pca",
                            new.reduction  = tolower(method),
                            group.by.vars  = unlist(int$group_by),
                            verbose        = FALSE),
                       extra)
                   )

                   merged <- Seurat::RunUMAP(
                     merged,
                     reduction      = tolower(method),
                     dims           = 1:pcs,
                     reduction.name = "umap",
                     verbose        = FALSE
                   )
                 }

                 # Join RNA layers — JoinLayers does not support SCTAssay
                 merged <- SeuratObject::JoinLayers(merged, assay = "RNA")
                 
                 if ("decontX" %in% names(merged@assays)) {
                   tryCatch(
                     merged <- SeuratObject::JoinLayers(merged, assay = "decontX"),
                     error = function(e) {
                       cli::cli_alert_warning(
                         "Could not join decontX layers: {conditionMessage(e)}"
                       )
                     }
                   )
                 }
                 
                 # Set default assay for downstream analyses
                 default_assay <- dplyr::case_when(
                   "SCT"     %in% names(merged@assays) ~ "SCT",
                   "decontX" %in% names(merged@assays) &&
                     isTRUE(cfg$decontX$use_decontX_counts) ~ "decontX",
                   TRUE ~ "RNA"
                 )
                 Seurat::DefaultAssay(merged) <- default_assay
                 cli::cli_alert_info("DefaultAssay set to: {default_assay}")

                 cli::cli_alert_success("Integration + UMAP done")
                 merged
               })
}
