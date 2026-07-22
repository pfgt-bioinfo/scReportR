# =============================================================================
# R/load.R
# 10x Genomics data loading.
# =============================================================================

#' Load all samples defined in params.yml
#'
#' Reads each `filtered_feature_bc_matrix` directory (or `.h5` file) and
#' creates an annotated Seurat v5 object per sample. Results are cached.
#'
#' @param cfg    Config list from [sc_config()].
#' @param force  If `TRUE`, ignore cache and reload from disk.
#'
#' @return A named list of Seurat objects (unfiltered).
#' @export
sc_load <- function(cfg, force = FALSE) {

  sc_cache_run("01_raw_objects.qs2",
               cache_dir = sc_cache_dir(cfg),
               force     = force,
               expr = {
                 meta <- sc_sample_meta(cfg)
                 cli::cli_h1("Loading {nrow(meta)} sample(s)")

                 obj_list <- lapply(seq_len(nrow(meta)), function(i) {
                   .load_one_sample(meta[i, ])
                 })
                 names(obj_list) <- meta$name

                 counts <- sapply(obj_list, ncol)
                 cli::cli_alert_success(
                   "Loaded: {paste(names(counts), counts, sep = ' = ', collapse = ', ')} cells"
                 )
                 obj_list
               })
}

# ── Internal ──────────────────────────────────────────────────────────────────

#' @keywords internal
.load_one_sample <- function(sample_info) {
  cli::cli_alert_info("Loading: {.strong {sample_info$name}}")
  
  path <- sample_info$path
  if (grepl("\\.h5$", path)) {
    mat <- Seurat::Read10X_h5(path)
  } else {
    mat <- Seurat::Read10X(data.dir = path)
  }
  
  if (is.list(mat)) mat <- mat[["Gene Expression"]]
  
  obj <- Seurat::CreateSeuratObject(
    counts       = mat,
    project      = sample_info$name,
    min.cells    = 3L,
    min.features = 10L
  )
  
  obj$sample_id <- sample_info$name
  
  # Everything else in the row is carried over as-is: the core descriptive
  # fields (library_type, species) plus whatever this project declared under
  # metadata: in params.yml. name and the paths are not metadata.
  fields <- setdiff(names(sample_info), c("name", "path", "raw_path"))
  
  clash <- intersect(fields, c("orig.ident", "nCount_RNA", "nFeature_RNA",
                               "sample_id", "seurat_clusters"))
  if (length(clash)) {
    cli::cli_abort(c(
      "Field name{?s} reserved by Seurat: {.field {clash}}.",
      "i" = "Rename {?it/them} in {.file params.yml}."
    ))
  }
  
  for (f in fields) obj[[f]] <- sample_info[[f]]
  
  obj <- scCustomize::Add_Cell_QC_Metrics(
    obj,
    species         = sample_info$species %||% "human",
    mito_name       = "percent.mt",
    ribo_name       = "percent.rb",
    mito_ribo_name  = "percent_mito_ribo",
    complexity_name = "log10_genes_per_umi",
    hemo_name       = "percent.hemo"
  )
  
  obj
}

#' Load empty droplets from a raw matrix (for decontX background)
#'
#' Reads `raw_feature_bc_matrix` and returns only barcodes absent from the
#' filtered matrix — these are true empty droplets used as ambient RNA
#' background for decontX.
#'
#' @param sample_info   One row of the sample metadata data frame.
#' @param filtered_mat  The filtered count matrix (dgCMatrix).
#'
#' @return A dgCMatrix of empty droplets, or `NULL` if `raw_path` is absent.
#' @export
sc_load_background <- function(sample_info, filtered_mat) {

  raw_path <- sample_info$raw_path %||% NA_character_

  if (is.na(raw_path) || !nzchar(raw_path)) {
    cli::cli_alert_warning(
      "No raw_path for {sample_info$name} — decontX will run without background"
    )
    return(NULL)
  }

  if (!file.exists(raw_path) && !dir.exists(raw_path)) {
    cli::cli_alert_warning(
      "raw_path not found: {.path {raw_path}} — running without background"
    )
    return(NULL)
  }

  cli::cli_alert_info("Loading background: {.path {raw_path}}")

  if (grepl("\\.h5$", raw_path)) {
    raw_mat <- Seurat::Read10X_h5(raw_path)
  } else {
    raw_mat <- Seurat::Read10X(data.dir = raw_path)
  }
  if (is.list(raw_mat)) raw_mat <- raw_mat[["Gene Expression"]]

  # Keep only empty droplets
  empty_bc <- setdiff(colnames(raw_mat), colnames(filtered_mat))

  if (length(empty_bc) == 0L) {
    cli::cli_alert_warning("No empty droplets found for {sample_info$name}")
    return(NULL)
  }

  common_genes <- intersect(rownames(raw_mat), rownames(filtered_mat))
  bg           <- raw_mat[common_genes, empty_bc, drop = FALSE]

  cli::cli_alert_success(
    "Background: {scales::comma(ncol(bg))} empty droplets"
  )
  bg
}
