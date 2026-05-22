# =============================================================================
# R/config.R
# Project configuration loader.
# =============================================================================

#' Load project configuration from params.yml
#'
#' Reads `params.yml`, creates output/cache directories, sets the random seed,
#' registers BiocParallel workers, and returns a named list of all parameters.
#' This should be the first call in `analysis.qmd`.
#'
#' @param params_file Path to the `params.yml` file. Defaults to `"params.yml"`
#'   in the current working directory.
#'
#' @return A named list of project parameters (the full parsed YAML).
#' @export
#'
#' @examples
#' \dontrun{
#' cfg <- sc_config("params.yml")
#' }
sc_config <- function(params_file = "params.yml") {

  if (!file.exists(params_file)) {
    cli::cli_abort(
      "params.yml not found at {.path {params_file}}.
       Run {.fn init_sc_report} to initialise the project."
    )
  }

  cfg <- yaml::read_yaml(params_file)

  # ── Create directories ──────────────────────────────────────────────────────
  dirs <- c(
    cfg$project$output_dir,
    cfg$project$cache_dir,
    file.path(cfg$project$output_dir, "figures")
  )
  invisible(lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE))

  # ── Random seed ─────────────────────────────────────────────────────────────
  set.seed(cfg$project$seed %||% 42L)

  # ── BiocParallel ────────────────────────────────────────────────────────────
  n_workers <- min(4L, max(1L, parallel::detectCores() - 1L))
  BiocParallel::register(BiocParallel::MulticoreParam(n_workers))

  # ── Seurat v5 ───────────────────────────────────────────────────────────────
  options(Seurat.object.assay.version = "v5")

  # ── Sample metadata as data frame ───────────────────────────────────────────
  cfg$.sample_meta <- dplyr::bind_rows(
    lapply(cfg$samples, function(s) {
      # fill optional fields with NA
      s$raw_path   <- s$raw_path   %||% NA_character_
      s$population <- s$population %||% NA_character_
      s$species    <- s$species    %||% "human"
      as.data.frame(s, stringsAsFactors = FALSE)
    })
  )

  # ── Colour palettes ───────────────────────────────────────────────────────────
  # All palettes can be overridden in params.yml under "colors:"
  # Defaults are auto-generated from the sample metadata.
  
  user_colors <- cfg$colors %||% list()
  
  # Tissues
  tissues <- unique(cfg$.sample_meta$tissue)
  tissue_defaults <- stats::setNames(
    RColorBrewer::brewer.pal(max(3L, length(tissues)), "Set1")[seq_along(tissues)],
    tissues
  )
  
  # Library types
  lib_types <- unique(cfg$.sample_meta$library_type)
  lib_defaults <- stats::setNames(
    RColorBrewer::brewer.pal(max(3L, length(lib_types)), "Dark2")[seq_along(lib_types)],
    lib_types
  )
  
  # Populations
  populations <- unique(cfg$.sample_meta$population)
  populations <- populations[!is.na(populations)]
  pop_defaults <- stats::setNames(
    RColorBrewer::brewer.pal(max(3L, length(populations)), "Set2")[seq_along(populations)],
    populations
  )
  
  # Samples
  sample_names <- cfg$.sample_meta$name
  sample_defaults <- stats::setNames(
    scales::hue_pal()(length(sample_names)),
    sample_names
  )
  
  .to_list <- function(x) if (is.null(x)) list() else as.list(x)
  
  cfg$.colors <- list(
    tissue       = unlist(utils::modifyList(as.list(tissue_defaults),
                                            .to_list(user_colors$tissue))),
    library_type = unlist(utils::modifyList(as.list(lib_defaults),
                                            .to_list(user_colors$library_type))),
    population   = unlist(utils::modifyList(as.list(pop_defaults),
                                            .to_list(user_colors$population))),
    samples      = unlist(utils::modifyList(as.list(sample_defaults),
                                            .to_list(user_colors$samples)))
  )

  cli::cli_alert_success(
    "Config loaded — project: {.strong {cfg$project$name}}, \\
     {nrow(cfg$.sample_meta)} sample(s)"
  )

  cfg
}

#' Return the cache directory from a config object
#' @param cfg Config list returned by [sc_config()].
#' @return Character path.
#' @export
sc_cache_dir <- function(cfg) cfg$project$cache_dir

#' Return the output directory from a config object
#' @param cfg Config list returned by [sc_config()].
#' @return Character path.
#' @export
sc_output_dir <- function(cfg) cfg$project$output_dir

#' Return sample metadata as a data frame
#' @param cfg Config list returned by [sc_config()].
#' @return A data frame with one row per sample.
#' @export
sc_sample_meta <- function(cfg) cfg$.sample_meta
