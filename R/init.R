# =============================================================================
# R/init.R
# Project initialisation — creates directory structure, params.yml, analysis.qmd
# =============================================================================

#' Initialise a new scRNAseq analysis project
#'
#' Creates the project directory with a pre-filled `params.yml`,
#' a ready-to-render `analysis.qmd`, and an RStudio project file.
#'
#' @param project_name  Short name for the project (used in filenames and titles).
#' @param path          Where to create the project directory. Defaults to the
#'   current working directory.
#' @param author        Author name for the report.
#' @param samples       Optional character vector of sample names to pre-fill
#'   in `params.yml`.
#' @param open          If `TRUE` and running in RStudio, open the project.
#'
#' @return Path to the created project directory (invisibly).
#' @export
#'
#' @examples
#' \dontrun{
#' init_sc_report(
#'   project_name = "Patient01",
#'   path         = "~/analyses",
#'   author       = "Jane Smith",
#'   samples      = c("CD45neg_T1", "CD45neg_T2", "CD45pos_T1")
#' )
#' }
init_sc_report <- function(project_name,
                            path    = ".",
                            author  = "Author",
                            samples = NULL,
                            open    = rlang::is_interactive()) {

  proj_dir <- file.path(path, project_name)

  if (dir.exists(proj_dir)) {
    cli::cli_abort(
      "Directory already exists: {.path {proj_dir}}.
       Choose a different name or path."
    )
  }

  # ── Create directories ──────────────────────────────────────────────────────
  dirs <- c(proj_dir,
            file.path(proj_dir, "data"),
            file.path(proj_dir, "cache"),
            file.path(proj_dir, "output"),
            file.path(proj_dir, "output", "figures"))
  invisible(lapply(dirs, dir.create, recursive = TRUE))

  # ── Copy and fill templates ─────────────────────────────────────────────────
  tmpl_dir <- system.file("templates", package = "scReportR")

  .fill_template(
    src  = file.path(tmpl_dir, "params.yml"),
    dest = file.path(proj_dir, "params.yml"),
    vars = list(
      PROJECT_NAME = project_name,
      AUTHOR       = author,
      SAMPLES_BLOCK = .build_samples_block(samples)
    )
  )

  .fill_template(
    src  = file.path(tmpl_dir, "analysis.qmd"),
    dest = file.path(proj_dir, "analysis.qmd"),
    vars = list(
      PROJECT_NAME = project_name,
      AUTHOR       = author
    )
  )
  
  # ── Copy styles directory ───────────────────────────────────────────────────
  styles_src  <- file.path(tmpl_dir, "styles")
  styles_dest <- file.path(proj_dir, "styles")
  
  if (dir.exists(styles_src)) {
    dir.create(styles_dest, showWarnings = FALSE)
    files <- list.files(styles_src, full.names = TRUE)
    file.copy(files, styles_dest)
    cli::cli_alert_success("Styles copied: {length(files)} file(s)")
  }

  # ── RStudio project file ────────────────────────────────────────────────────
  rproj_content <- paste0(
    "Version: 1.0\n\n",
    "RestoreWorkspace: No\n",
    "SaveWorkspace: No\n",
    "AlwaysSaveHistory: No\n\n",
    "EnableCodeIndexing: Yes\n",
    "UseSpacesForTab: Yes\n",
    "NumSpacesForTab: 2\n",
    "Encoding: UTF-8\n\n",
    "RnwWeave: knitr\n",
    "LaTeX: pdfLaTeX\n"
  )
  writeLines(rproj_content,
             file.path(proj_dir, paste0(project_name, ".Rproj")))

  # ── Summary ─────────────────────────────────────────────────────────────────
  cli::cli_h1("Project initialised: {project_name}")
  cli::cli_ul(c(
    "Directory  : {.path {proj_dir}}",
    "Edit       : {.path {file.path(proj_dir, 'params.yml')}}",
    "Render     : {.code quarto::quarto_render('analysis.qmd')}"
  ))

  if (open && requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    rstudioapi::openProject(proj_dir, newSession = TRUE)
  }

  invisible(proj_dir)
}

# ── Internal helpers ──────────────────────────────────────────────────────────

#' @keywords internal
.fill_template <- function(src, dest, vars) {
  content <- paste(readLines(src, warn = FALSE), collapse = "\n")
  for (nm in names(vars)) {
    content <- gsub(paste0("\\{\\{", nm, "\\}\\}"), vars[[nm]], content)
  }
  writeLines(content, dest)
}

#' @keywords internal
.build_samples_block <- function(samples) {
  if (is.null(samples) || length(samples) == 0L) {
    return(paste0(
      "  - name: \"sample_01\"\n",
      "    path: \"data/sample_01/filtered_feature_bc_matrix\"\n",
      "    raw_path: \"data/sample_01/raw_feature_bc_matrix\"\n",
      "    library_type: \"OCM\"   # OCM | 3pv4 | 5p | multiome\n",
      "    tissue: \"T1\"\n",
      "    population: \"CD45neg\"\n"
    ))
  }

  paste(sapply(seq_along(samples), function(i) {
    nm <- samples[i]
    paste0(
      "  - name: \"", nm, "\"\n",
      "    path: \"data/", nm, "/filtered_feature_bc_matrix\"\n",
      "    raw_path: \"data/", nm, "/raw_feature_bc_matrix\"\n",
      "    library_type: \"\"   # OCM | 3pv4 | 5p | multiome\n",
      "    tissue: \"\"\n",
      "    population: \"\"\n"
    )
  }), collapse = "\n")
}
