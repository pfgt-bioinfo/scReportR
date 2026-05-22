# =============================================================================
# R/qc.R
# QC metric computation and threshold resolution.
# No plotting — all visualisation lives in analysis.qmd.
# =============================================================================

#' Get QC thresholds for a sample (defaults + per-sample overrides)
#'
#' @param sample_name  Sample identifier.
#' @param cfg          Config list from [sc_config()].
#'
#' @return Named list of QC thresholds.
#' @export
sc_qc_thresholds <- function(sample_name, cfg) {
  thr      <- cfg$qc$defaults
  override <- cfg$qc$per_sample[[sample_name]]
  if (!is.null(override)) {
    for (nm in names(override)) thr[[nm]] <- override[[nm]]
  }
  thr
}

#' Get thresholds for all samples
#'
#' @param cfg  Config list from [sc_config()].
#' @return Named list of threshold lists, one per sample.
#' @export
sc_qc_all_thresholds <- function(cfg) {
  sample_names <- sc_sample_meta(cfg)$name
  lapply(sample_names, sc_qc_thresholds, cfg = cfg) |>
    stats::setNames(sample_names)
}

#' Compute a pre-filtering QC summary table
#'
#' Optional metadata columns (tissue, population) are included only when
#' present in the object — the function never fails on missing columns.
#'
#' @param obj_list  Named list of Seurat objects.
#' @return A tibble with one row per sample.
#' @export
sc_qc_summary <- function(obj_list) {
  dplyr::bind_rows(lapply(names(obj_list), function(sname) {
    obj  <- obj_list[[sname]]
    meta <- obj@meta.data

    # Core metrics — always present
    row <- tibble::tibble(
      sample        = sname,
      library_type  = .meta_unique(meta, "library_type"),
      n_cells       = ncol(obj),
      median_genes  = round(stats::median(meta$nFeature_RNA)),
      median_umi    = round(stats::median(meta$nCount_RNA)),
      median_mt_pct = round(stats::median(meta$percent.mt), 1),
      pct_mt20      = round(100 * mean(meta$percent.mt > 20), 1)
    )

    # Optional metadata columns — included only when available
    for (col in c("tissue", "population", "species")) {
      if (col %in% colnames(meta)) {
        row[[col]] <- .meta_unique(meta, col)
      }
    }
    row
  }))
}

#' Build a data frame of thresholds (for plotting)
#'
#' @param thresholds_list  Output of [sc_qc_all_thresholds()].
#' @param cfg              Config list.
#' @return A tibble with one row per sample, including override flag.
#' @export
sc_qc_thresholds_df <- function(thresholds_list, cfg) {
  dplyr::bind_rows(lapply(names(thresholds_list), function(nm) {
    thr <- thresholds_list[[nm]]
    tibble::tibble(
      sample       = nm,
      has_override = nm %in% names(cfg$qc$per_sample %||% list()),
      !!!thr
    )
  }))
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' Return a single value from a metadata column, collapsing multiple values
#' @keywords internal
.meta_unique <- function(meta, col) {
  if (!col %in% colnames(meta)) return(NA_character_)
  vals <- unique(stats::na.omit(meta[[col]]))
  if (length(vals) == 0L) return(NA_character_)
  if (length(vals) == 1L) return(as.character(vals))
  paste(vals, collapse = "/")
}

#' Resolve a colour palette for a given metadata column
#'
#' Returns the palette from cfg if it matches the column name,
#' otherwise returns NULL (letting ggplot pick defaults).
#'
#' @param col  Metadata column name.
#' @param cfg  Config list from [sc_config()].
#' @keywords internal
.resolve_palette <- function(col, cfg) {
  if (is.null(cfg)) return(NULL)
  switch(col,
    sample_id    = cfg$.colors$samples,
    library_type = cfg$.colors$library_type,
    tissue       = cfg$.colors$tissue,
    population   = cfg$.colors$population,
    NULL
  )
}

#' Pick the best available categorical metadata column for fill/colour
#'
#' Tries columns in order; returns the first one present.
#'
#' @param meta      Metadata data frame.
#' @param preferred Character vector of preferred column names, in priority order.
#' @keywords internal
.best_meta_col <- function(meta, preferred = c("library_type", "tissue",
                                                "population", "sample_id")) {
  found <- preferred[preferred %in% colnames(meta)]
  if (length(found) == 0L) return(NULL)
  found[[1L]]
}
