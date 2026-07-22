# =============================================================================
# R/cluster.R
# Clustering (Leiden/Louvain) + differential markers.
# =============================================================================

#' Cluster cells at multiple resolutions
#'
#' Runs FindNeighbors then FindClusters at all resolutions defined in
#' `params.yml`. Sets Idents to the default resolution.
#'
#' @param merged  Integrated Seurat object from [sc_integrate()].
#' @param cfg     Config list from [sc_config()].
#' @param force   If `TRUE`, ignore cache.
#'
#' @return Clustered Seurat object.
#' @export
sc_cluster <- function(merged, cfg, force = FALSE) {

  sc_cache_run("04a_clustered.qs2",
               cache_dir = sc_cache_dir(cfg),
               force     = force,
               expr = {
                 clust  <- cfg$clustering
                 method <- cfg$integration$method %||% "Harmony"
                 pcs    <- cfg$normalization$pcs_use %||% 30L

                 reduction <- switch(method,
                   Harmony         = "harmony",
                   RPCAIntegration = "rpcaintegration",
                   CCAIntegration  = "ccaintegration",
                   none            = "pca"
                 )

                 cli::cli_h1("Clustering (algo={clust$algorithm %||% 4})")

                 merged <- Seurat::FindNeighbors(
                   merged, reduction = reduction,
                   dims = 1:pcs, verbose = FALSE
                 )

                 for (res in clust$resolutions) {
                   col <- paste0("clusters_res", res)
                   cli::cli_alert_info("Resolution: {res}")
                   merged <- Seurat::FindClusters(
                     merged,
                     resolution   = res,
                     algorithm    = clust$algorithm %||% 4L,
                     random.seed = 1, 
                     cluster.name = col,
                     verbose      = FALSE
                   )
                 }

                 default_col <- paste0("clusters_res",
                                       clust$default_resolution %||% 0.6)
                 Seurat::Idents(merged) <- default_col

                 cli::cli_alert_success(
                   "Clustering done — {length(levels(Seurat::Idents(merged)))} \\
                    clusters at res={clust$default_resolution %||% 0.6}"
                 )
                 merged
               })
}

#' Find differential markers for all clusters
#'
#' Automatically uses SCT assay if available and calls PrepSCTFindMarkers
#' when needed. Uses parameters from `params.yml`.
#'
#' @param merged  Clustered Seurat object from [sc_cluster()].
#' @param cfg     Config list from [sc_config()].
#' @param n_top   Number of top markers per cluster to return.
#' @param force   If `TRUE`, ignore cache.
#'
#' @return A data frame of markers (output of FindAllMarkers).
#' @export
sc_find_markers <- function(merged, cfg, n_top = 5L, force = FALSE) {

  sc_cache_run("04b_markers.qs2",
               cache_dir = sc_cache_dir(cfg),
               force     = force,
               expr = {
                 assay       <- if ("SCT" %in% names(merged@assays)) "SCT" else "RNA"
                 cluster_col <- paste0("clusters_res",
                                       cfg$clustering$default_resolution %||% 0.6)
                 min_pct     <- cfg$clustering$findmarkers_min_pct %||% 0.1
                 logfc       <- cfg$clustering$findmarkers_logfc   %||% 0.25

                 Seurat::Idents(merged)       <- cluster_col
                 Seurat::DefaultAssay(merged) <- assay

                 if (assay == "SCT") {
                   cli::cli_alert_info("Running PrepSCTFindMarkers...")
                   future::plan("sequential")
                   BiocParallel::register(BiocParallel::SerialParam())
                   options(future.globals.maxSize = 8 * 1024^3)
                   merged <- Seurat::PrepSCTFindMarkers(merged, verbose = FALSE)
                   future::plan("multisession")
                   BiocParallel::register(
                     BiocParallel::MulticoreParam(
                       min(4L, parallel::detectCores() - 1L)
                     )
                   )
                   options(future.globals.maxSize = NULL)
                 }

                 cli::cli_h1("FindAllMarkers (assay={assay}, min.pct={min_pct})")

                 markers <- Seurat::FindAllMarkers(
                   merged,
                   only.pos        = TRUE,
                   min.pct         = min_pct,
                   logfc.threshold = logfc,
                   layer           = "data",
                   test.use        = "wilcox",
                   verbose         = FALSE
                 )

                 if (nrow(markers) == 0L || !"cluster" %in% colnames(markers)) {
                   cli::cli_alert_warning("No significant markers found")
                 } else {
                   cli::cli_alert_success(
                     "{nrow(markers)} markers across \\
                      {length(unique(markers$cluster))} clusters"
                   )
                 }
                 markers
               })
}

#' Compute module scores for top markers of each cluster
#'
#' Uses AddModuleScore to score each cluster's top markers on all cells,
#' then stores results in metadata as "score_cluster_X" columns.
#'
#' @param merged      Clustered Seurat object.
#' @param markers_df  Output of [sc_find_markers()].
#' @param cfg         Config list from [sc_config()].
#' @param n_top       Number of top markers per cluster to use.
#' @param force       If TRUE, ignore cache.
#'
#' @return Seurat object with module score columns added to metadata.
#' @export
sc_module_scores <- function(merged, markers_df, cfg, n_top = 20L, force = FALSE) {
  
  sc_cache_run("04c_module_scores.qs2",
               cache_dir = sc_cache_dir(cfg),
               force     = force,
               expr = {
                 if (is.null(markers_df) || nrow(markers_df) == 0L) {
                   cli::cli_alert_warning("No markers available — skipping module scores")
                   return(merged)
                 }
                 
                 clusters <- sort(unique(markers_df$cluster))
                 cli::cli_h1("Computing module scores for {length(clusters)} clusters")
                 
                 for (clust in clusters) {
                   genes <- markers_df |>
                     dplyr::filter(cluster == clust) |>
                     dplyr::slice_max(avg_log2FC, n = n_top, with_ties = FALSE) |>
                     dplyr::pull(gene)
                   
                   if (length(genes) < 3L) {
                     cli::cli_alert_warning("Cluster {clust}: fewer than 3 markers, skipping")
                     next
                   }
                   
                   col_name <- paste0("score_cluster_", clust)
                   cli::cli_alert_info("Cluster {clust}: {length(genes)} genes")
                   
                   merged <- Seurat::AddModuleScore(
                     merged,
                     features = list(genes),
                     name     = col_name,
                     verbose  = FALSE
                   )
                   
                   # AddModuleScore appends "1" to the name — rename to keep clean
                   merged@meta.data[[col_name]] <- merged@meta.data[[paste0(col_name, "1")]]
                   merged@meta.data[[paste0(col_name, "1")]] <- NULL
                 }
                 
                 cli::cli_alert_success("Module scores computed")
                 merged
               })
}