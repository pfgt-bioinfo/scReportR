# =============================================================================
# R/cache.R
# qs2-based cache system with force-recompute support.
# =============================================================================

#' Load from cache or compute and save
#'
#' If a `.qs2` cache file exists, loads and returns it. Otherwise evaluates
#' `expr`, saves the result, and returns it. Designed to wrap expensive
#' computation steps (normalisation, clustering, etc.).
#'
#' @param cache_file  File name (without path) of the `.qs2` cache file.
#' @param expr        Expression to evaluate if cache is missing.
#' @param cache_dir   Directory where cache files are stored.
#' @param force       If `TRUE`, recompute even if cache exists.
#'
#' @return The cached or freshly computed result.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- sc_cache_run("my_step.qs2", expr = {
#'   Sys.sleep(2)
#'   42
#' }, cache_dir = "cache")
#' }
sc_cache_run <- function(cache_file, expr, cache_dir, force = FALSE) {

  path <- file.path(cache_dir, cache_file)

  if (file.exists(path) && !isTRUE(force)) {
    cli::cli_alert_info("Cache hit  : {.file {cache_file}}")
    return(qs2::qs_read(path))
  }

  cli::cli_alert_info("Computing  : {.file {cache_file}}")
  result <- eval(expr, envir = parent.frame())
  qs2::qs_save(result, path)
  cli::cli_alert_success("Saved      : {.file {cache_file}}")
  invisible(result)
}

#' List all cache files for a project
#'
#' @param cache_dir  Cache directory path.
#' @return Character vector of `.qs2` file paths.
#' @export
sc_cache_list <- function(cache_dir) {
  list.files(cache_dir, pattern = "\\.qs2$", full.names = TRUE)
}

#' Delete one or all cache files
#'
#' @param cache_dir   Cache directory path.
#' @param cache_file  Optional. Specific file name to delete. If `NULL`,
#'   all cache files are deleted.
#' @export
sc_cache_clear <- function(cache_dir, cache_file = NULL) {
  if (is.null(cache_file)) {
    files <- sc_cache_list(cache_dir)
    file.remove(files)
    cli::cli_alert_success("Cleared {length(files)} cache file(s).")
  } else {
    path <- file.path(cache_dir, cache_file)
    if (file.exists(path)) {
      file.remove(path)
      cli::cli_alert_success("Deleted: {.file {cache_file}}")
    } else {
      cli::cli_alert_warning("Cache file not found: {.file {cache_file}}")
    }
  }
  invisible(NULL)
}

#' Keep only essential cache files and remove intermediates
#'
#' After the analysis is finalised, removes all intermediate cache files
#' to save disk space. Only the raw objects and the final annotated object
#' are kept — the full analysis can be reproduced from these two files.
#'
#' @param cfg      Config list from [sc_config()].
#' @param dry_run  If `TRUE`, print what would be deleted without deleting.
#'   Defaults to `TRUE` for safety.
#' @export
sc_cache_finalize <- function(cfg, dry_run = TRUE) {

  keep  <- c("01_raw_objects.qs2", "final_annotated.qs2")
  all_f <- sc_cache_list(sc_cache_dir(cfg))
  to_keep   <- all_f[basename(all_f) %in% keep]
  to_delete <- all_f[!basename(all_f) %in% keep]

  cli::cli_h1("Cache finalisation")

  cli::cli_alert_success("Files to keep ({length(to_keep)}):")
  lapply(to_keep, function(f) cli::cli_alert_success("  {basename(f)}"))

  cli::cli_alert_warning("Files to remove ({length(to_delete)}):")
  lapply(to_delete, function(f) cli::cli_alert_warning("  {basename(f)}"))

  if (dry_run) {
    cli::cli_alert_info(
      "Dry run — nothing deleted. Run with {.code dry_run = FALSE} to confirm."
    )
    return(invisible(NULL))
  }

  file.remove(to_delete)
  cli::cli_alert_success(
    "Done — {length(to_delete)} intermediate file(s) removed."
  )
  invisible(NULL)
}
