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
      PROJECT_NAME   = project_name,
      AUTHOR         = author,
      DEFAULTS_BLOCK = .build_defaults_block(), 
      SAMPLES_BLOCK  = .build_samples_block(samples),
      METADATA_BLOCK = .build_metadata_block()
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
    # Recursive: styles/fonts/ holds the .woff2 that pfgt_theme.scss declares
    # via @font-face, and embed-resources inlines them at render time — a
    # missing font is a hard Quarto error, not a fallback.
    files <- list.files(styles_src, recursive = TRUE)   # relative paths, files only
    
    subdirs <- setdiff(unique(dirname(files)), ".")
    for (d in subdirs) {
      dir.create(file.path(styles_dest, d), recursive = TRUE, showWarnings = FALSE)
    }
    dir.create(styles_dest, showWarnings = FALSE)
    
    ok <- file.copy(
      file.path(styles_src, files),
      file.path(styles_dest, files),
      overwrite = TRUE
    )
    
    # file.copy() never errors — it returns FALSE. Say so, or the next missing
    # asset will surface as a Quarto error three steps later.
    if (!all(ok)) {
      cli::cli_abort(c(
        "Could not copy {sum(!ok)} style asset{?s}:",
        "x" = "{.file {files[!ok]}}"
      ))
    }
    cli::cli_alert_success("Styles copied: {length(files)} file{?s}")
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
    content <- gsub(paste0("{{", nm, "}}"), vars[[nm]], content, fixed = TRUE)
  }
  writeLines(content, dest)
}

#' Build the `defaults:` block of a starter params.yml
#'
#' Fields shared by every sample. sc_config() merges these into each sample,
#' where a per-sample value takes precedence.
#'
#' @keywords internal
.build_defaults_block <- function() {
  paste0(
    '  library_type: "OCM"   # OCM | 3pv4 | 5p | multiome | Flex\n',
    '  species: "human"\n'
  )
}

#' Build the `samples:` block of a starter params.yml
#'
#' @param samples Character vector of sample names, or `NULL` for a single
#'   placeholder entry.
#' @param metadata_fields Character vector of descriptive field names to stub
#'   as flat keys on each sample. Placeholders meant to be renamed per project —
#'   the package never refers to them by name. `NULL` omits them.
#'
#' @keywords internal
.build_samples_block <- function(samples,
                                 metadata_fields = c("condition", "timepoint")) {
  
  md <- if (length(metadata_fields)) {
    paste0(
      "    # Descriptive fields — flat keys, rename to suit this project.\n",
      "    # tumour: tissue, population | cell line: treatment, timepoint, coculture\n",
      "    # Quote every value: unquoted yes/no/on/off are read as booleans.\n",
      paste0('    ', metadata_fields, ': ""\n', collapse = "")
    )
  } else {
    ""
  }
  
  one <- function(nm) paste0(
    "  - name: \"", nm, "\"\n",
    "    path: \"data/", nm, "/filtered_feature_bc_matrix\"\n",
    "    raw_path: \"data/", nm, "/raw_feature_bc_matrix\"\n",
    md
  )
  
  if (is.null(samples) || length(samples) == 0L) return(one("sample_01"))
  paste(vapply(samples, one, character(1)), collapse = "\n")
}

#' Build the `metadata_labels:` block of a starter params.yml
#'
#' Emits one commented `field: "Field"` line per metadata field. Kept commented
#' so the block stays inert: labels are optional, and an active block naming the
#' placeholder fields would trip sc_config()'s stray-label warning the moment a
#' user renames their metadata: fields without updating the labels.
#'
#' @param metadata_fields Character vector, the same fields stubbed under each
#'   sample's `metadata:` by [.build_samples_block()]. `NULL` yields a single
#'   placeholder comment.
#'
#' @keywords internal
.build_metadata_block <- function(metadata_fields = c("condition", "timepoint")) {
  if (!length(metadata_fields)) return("  # condition: \"Condition\"\n")
  
  # "condition" -> "Condition": a readable default the user can edit or delete.
  labels <- paste0(
    toupper(substring(metadata_fields, 1L, 1L)),
    substring(metadata_fields, 2L)
  )
  paste0('  # ', metadata_fields, ': "', labels, '"\n', collapse = "")
}
