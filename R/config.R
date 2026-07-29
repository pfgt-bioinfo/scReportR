# =============================================================================
# R/config.R
# Project configuration loader.
# =============================================================================

# ── Internal constants ───────────────────────────────────────────────────────

#' Core sample fields recognised by `sc_config()`
#'
#' Fields that every project shares and that the code is allowed to know about
#' by name. Any other flat key on a sample is treated as a free-form
#' descriptive field instead: discovered at load time, added to
#' `cfg$.sample_meta`, and given an automatic colour palette.
#'
#' @noRd
.SC_CORE_FIELDS <- c("name", "path", "raw_path", "library_type", "species")

#' Legacy top-level fields promoted into `metadata:` for backward compatibility
#'
#' Older `params.yml` files declared these directly on the sample. They are
#' moved into `metadata:` with a warning so that existing projects keep working.
#'
#' @noRd
.SC_LEGACY_FIELDS <- c("tissue", "population")

#' Default ColorBrewer palettes, cycled across discovered metadata fields
#'
#' Assigned in the order fields appear in `cfg$.meta_fields`, and recycled if a
#' project declares more fields than there are palettes. Override per field via
#' `palettes:` in `params.yml`.
#'
#' @noRd
.SC_DEFAULT_PALETTES <- c("Set1", "Dark2", "Set2", "Paired", "Accent", "Set3")


# ── Internal helpers ─────────────────────────────────────────────────────────

#' Coerce a value to a plain list
#'
#' Convenience wrapper so that absent optional YAML blocks (`colors:`,
#' `palettes:`, `metadata_labels:`) can be passed straight to
#' [utils::modifyList()], which requires a list on both sides.
#'
#' @param x Any object, typically a list read from YAML, or `NULL`.
#'
#' @return `x` as a list; `list()` when `x` is `NULL`.
#'
#' @noRd
.sc_to_list <- function(x) if (is.null(x)) list() else as.list(x)

#' Build a named colour vector for a set of levels
#'
#' Unlike a bare [RColorBrewer::brewer.pal()] call, this never silently returns
#' fewer colours than there are levels: past the palette's maximum (8 for most
#' qualitative palettes, 9 for `"Set1"`, 12 for `"Set3"`) it interpolates with
#' [grDevices::colorRampPalette()] rather than dropping the extra levels.
#'
#' @param values Vector of levels. Duplicates and `NA` are dropped; the
#'   remaining values are used as the names of the returned vector, in order of
#'   first appearance.
#' @param palette Name of a ColorBrewer palette, as listed in
#'   [RColorBrewer::brewer.pal.info].
#'
#' @return A named character vector of hex colours, one per level, or `NULL`
#'   when `values` holds no non-missing level.
#'
#' @noRd
.sc_palette_for <- function(values, palette = "Set1") {
  values <- unique(values[!is.na(values)])
  n <- length(values)
  if (n == 0L) return(NULL)
  
  if (!palette %in% rownames(RColorBrewer::brewer.pal.info)) {
    cli::cli_abort("Unknown ColorBrewer palette: {.val {palette}}.")
  }
  max_n <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
  
  cols <- if (n <= max_n) {
    RColorBrewer::brewer.pal(max(3L, n), palette)[seq_len(n)]
  } else {
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(max_n, palette))(n)
  }
  stats::setNames(cols, values)
}

#' Flatten one `cfg$samples` entry into a data frame
#'
#' Merges the project-level `defaults:` block into the sample (a per-sample
#' value overrides the default), then flattens all keys. Descriptive fields are
#' flat, alongside the core fields — anything that isn't core is discovered as
#' metadata downstream. Fields whose length is neither 1 nor `length(name)`
#' raise an error naming the offender, rather than being silently recycled.
#'
#' @param s        One entry of `cfg$samples`, as read from YAML.
#' @param defaults The parsed `defaults:` block (list), applied to every sample.
#'
#' @return A data frame with one row per element of `s$name`.
#'
#' @noRd
.sc_sample_row <- function(s, defaults = list()) {
  
  # Project-level defaults fill fields the sample omits; a per-sample value
  # takes precedence (modifyList overwrites `defaults` with `s`).
  s <- utils::modifyList(defaults, s)
  
  # Optional core fields still get a hard default if neither sample nor
  # defaults: supplied them.
  s$raw_path <- s$raw_path %||% NA_character_
  s$species  <- s$species  %||% "human"
  
  if (is.null(s$name)) {
    cli::cli_abort("A sample entry has no {.field name}.")
  }
  
  flat <- lapply(s, function(v) unlist(v, use.names = FALSE))
  
  n <- length(flat$name)
  bad <- names(flat)[!lengths(flat) %in% c(1L, n)]
  if (length(bad)) {
    cli::cli_abort(c(
      "Length mismatch in sample {.val {flat$name}}:",
      "x" = "{.field {bad}} — expected 1 or {n} value{?s}."
    ))
  }
  
  flat <- lapply(flat, function(v) if (length(v) == 1L) rep(v, n) else v)
  as.data.frame(flat, stringsAsFactors = FALSE)
}

#' Safely walk a nested path in the parsed config
#'
#' Returns `NULL` for any missing key, and for a path that runs into a
#' non-list — unlike `[[`, which errors on an atomic vector.
#'
#' @noRd
.sc_get <- function(cfg, path) {
  x <- cfg
  for (n in path) {
    if (!is.list(x)) return(NULL)
    x <- x[[n]]
    if (is.null(x)) return(NULL)
  }
  x
}


# ── Main ─────────────────────────────────────────────────────────────────────

#' Load project configuration from params.yml
#'
#' Reads `params.yml`, creates output/cache directories, sets the random seed,
#' registers BiocParallel workers, and returns a named list of all parameters.
#' This should be the first call in `analysis.qmd`.
#'
#' Descriptive sample fields are declared under each sample's `metadata:` block
#' and are **not** hard-coded: whatever keys you use there are discovered at load
#' time, flattened into `cfg$.sample_meta`, listed in `cfg$.meta_fields`, and
#' given a colour palette in `cfg$.colors`. This lets the same code serve a
#' tumour project (`tissue`, `population`) and a cell-line project
#' (`treatment`, `timepoint`, `coculture`) without modification.
#'
#' @param params_file Path to the `params.yml` file. Defaults to `"params.yml"`
#'   in the current working directory.
#'
#' @return The full parsed YAML as a named list, plus:
#'   \describe{
#'     \item{`.sample_meta`}{Data frame, one row per sample: the core fields
#'       plus one column per metadata field.}
#'     \item{`.meta_fields`}{Character vector of the metadata fields discovered
#'       in `params.yml`. Iterate over this rather than naming fields.}
#'     \item{`.labels`}{Named list of display labels, one per metadata field.
#'       Defaults to the field name; override via `metadata_labels:`.}
#'     \item{`.colors`}{Named list of colour vectors: one per metadata field,
#'       plus `library_type` and `samples`. Palettes are chosen automatically;
#'       override the whole palette via `palettes:` or individual levels via
#'       `colors:`.}
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' cfg <- sc_config("params.yml")
#'
#' # Iterate over whatever fields this project happens to declare
#' for (f in cfg$.meta_fields) {
#'   print(scCustomize::DimPlot_scCustom(
#'     merged, group.by = f, colors_use = cfg$.colors[[f]]
#'   ) + ggplot2::ggtitle(cfg$.labels[[f]]))
#' }
#' }
sc_config <- function(params_file = "params.yml") {
  
  if (!file.exists(params_file)) {
    cli::cli_abort(
      "params.yml not found at {.path {params_file}}.
       Run {.fn init_sc_report} to initialise the project."
    )
  }
  
  cfg <- yaml::read_yaml(params_file)
  
  if (!length(cfg$samples)) {
    cli::cli_abort("No {.field samples} declared in {.path {params_file}}.")
  }
  
  # ── Config validation ───────────────────────────────────────────────────────
  # YAML validates nothing: a mistyped key is read back faithfully and silently
  # (rresolutions, fale, ...), so the failure surfaces far from its cause.
  # Below: keys we cannot run without, then keys we don't recognise.
  required <- list(
    c("project", "name"),
    c("normalization", "method"),
    c("integration", "method"),
    c("clustering", "resolutions"), 
    c("clustering", "default_resolution")
  )
  
  absent <- Filter(function(p) is.null(.sc_get(cfg, p)), required)
  
  if (length(absent)) {
    cli::cli_abort(c(
      "Missing required key{?s} in {.file {params_file}}:",
      "x" = "{.field {vapply(absent, paste, character(1), collapse = ' > ')}}"
    ))
  }
  
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
  # Slurm allocates a subset of the node's cores; detectCores() reports the
  # whole machine. Honour the allocation when there is one, and leave a core
  # free only when running locally.
  n_avail <- suppressWarnings(
    as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))
  )
  if (is.na(n_avail)) {
    n_avail <- parallel::detectCores(logical = FALSE)
    if (is.na(n_avail)) n_avail <- 1L   # detectCores() can return NA
    n_avail <- n_avail - 1L             # local: keep the machine usable
  }
  n_workers <- min(cfg$project$max_workers %||% 4L, max(1L, n_avail))
  
  BiocParallel::register(BiocParallel::MulticoreParam(
    n_workers, RNGseed = cfg$project$seed %||% 42L
  ))
  
  # ── Seurat v5 ───────────────────────────────────────────────────────────────
  options(Seurat.object.assay.version = "v5")
  
  # ── Sample metadata as data frame ───────────────────────────────────────────
  cfg$.sample_meta <- dplyr::bind_rows(
    lapply(cfg$samples, .sc_sample_row, defaults = .sc_to_list(cfg$defaults))
  )
  
  if (anyDuplicated(cfg$.sample_meta$name)) {
    dup <- unique(cfg$.sample_meta$name[duplicated(cfg$.sample_meta$name)])
    cli::cli_abort("Duplicated sample name{?s}: {.val {dup}}.")
  }
  
  # ── Metadata fields ─────────────────────────────────────────────────────────
  # Anything that isn't a core field is, by definition, this project's metadata.
  cfg$.meta_fields <- setdiff(names(cfg$.sample_meta), .SC_CORE_FIELDS)
  
  # bind_rows() fills gaps with NA when samples declare different fields.
  ragged <- vapply(
    cfg$.meta_fields,
    function(f) anyNA(cfg$.sample_meta[[f]]),
    logical(1)
  )
  if (any(ragged)) {
    cli::cli_warn(c(
      "Metadata field{?s} missing for some samples:
       {.field {cfg$.meta_fields[ragged]}}.",
      "i" = "Those samples will carry {.code NA} for these fields."
    ))
  }
  
  # YAML 1.1 reads unquoted yes/no/on/off/y/n as booleans — almost never what a
  # descriptive field means. Flag it rather than let it surface as a legend
  # reading "FALSE".
  logical_fields <- cfg$.meta_fields[
    vapply(cfg$.sample_meta[cfg$.meta_fields], is.logical, logical(1))
  ]
  if (length(logical_fields)) {
    cli::cli_warn(c(
      "Metadata field{?s} parsed as logical: {.field {logical_fields}}.",
      "i" = "YAML reads unquoted {.val yes}/{.val no}/{.val on}/{.val off}
             as booleans.",
      ">" = "Quote the values in {.file {params_file}} to keep them as text."
    ))
  }
  
  # ── A run: flag that isn't logical is almost always a typo (fale, ture, yes).
  for (p in list(c("qc", "run_prefilter_umap"),
                 c("decontX", "run"),
                 c("annotation", "singler", "run"),
                 c("annotation", "panhuman_azimuth", "run"))) {
    v <- .sc_get(cfg, p)
    if (!is.null(v) && !is.logical(v)) {
      cli::cli_warn(c(
        "{.field {paste(p, collapse = '$')}} is {.cls {class(v)}}, not logical: {.val {v}}.",
        ">" = "Expected {.val true} or {.val false} in {.file {params_file}}."
      ))
    }
  }
  
  # ── Display labels ──────────────────────────────────────────────────────────
  # Default to the field name; metadata_labels: overrides any subset of them.
  user_labels <- .sc_to_list(cfg$metadata_labels)
  
  displayable <- c("sample_id", "library_type", "species", cfg$.meta_fields)
  
  stray_labels <- setdiff(names(user_labels), displayable)
  if (length(stray_labels)) {
    cli::cli_warn(c(
      "Labels for unknown field{?s}: {.field {stray_labels}}.",
      "i" = "Declared metadata fields: {.field {cfg$.meta_fields}}."
    ))
  }
  
  cfg$.labels <- utils::modifyList(
    stats::setNames(as.list(displayable), displayable),
    user_labels
  )
  
  # ── Colour palettes ─────────────────────────────────────────────────────────
  # Auto-generated from the sample metadata. Two override levels, both optional:
  #   palettes:  <field>: <ColorBrewer palette name>   — swap the whole palette
  #   colors:    <field>: {<level>: <hex>}             — pin individual levels
  user_colors <- .sc_to_list(cfg$colors)
  user_pals   <- .sc_to_list(cfg$palettes)
  
  cfg$.colors <- list()
  
  for (i in seq_along(cfg$.meta_fields)) {
    f <- cfg$.meta_fields[[i]]
    pal_name <- user_pals[[f]] %||%
      .SC_DEFAULT_PALETTES[[((i - 1L) %% length(.SC_DEFAULT_PALETTES)) + 1L]]
    
    defaults <- .sc_palette_for(cfg$.sample_meta[[f]], pal_name)
    if (is.null(defaults)) next
    
    cfg$.colors[[f]] <- unlist(utils::modifyList(
      as.list(defaults), .sc_to_list(user_colors[[f]])
    ))
  }
  
  # colors:/palettes: naming a field that doesn't exist is a silent no-op —
  # exactly what happens when metadata: fields get renamed but the colour
  # blocks don't follow.
  stray_colors <- setdiff(names(user_colors),
                          c(cfg$.meta_fields, "library_type", "samples"))
  stray_pals   <- setdiff(names(user_pals), cfg$.meta_fields)
  
  if (length(stray_colors) || length(stray_pals)) {
    cli::cli_warn(c(
      "Colour overrides for unknown field{?s}: {.field {unique(c(stray_colors, stray_pals))}}.",
      "i" = "Declared metadata fields: {.field {cfg$.meta_fields}}.",
      ">" = "{?This override is/These overrides are} ignored."
    ))
  }
  
  # library_type and samples are core fields, so they are always available.
  cfg$.colors$library_type <- unlist(utils::modifyList(
    as.list(.sc_palette_for(cfg$.sample_meta$library_type, "Dark2")),
    .sc_to_list(user_colors$library_type)
  ))
  
  cfg$.colors$samples <- unlist(utils::modifyList(
    as.list(stats::setNames(
      scales::hue_pal()(nrow(cfg$.sample_meta)),
      cfg$.sample_meta$name
    )),
    .sc_to_list(user_colors$samples)
  ))
  
  cli::cli_alert_success(
    "Config loaded — project: {.strong {cfg$project$name}}, \\
     {nrow(cfg$.sample_meta)} sample{?s}"
  )
  if (length(cfg$.meta_fields)) {
    cli::cli_alert_info("Metadata fields: {.field {cfg$.meta_fields}}")
  }
  
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
