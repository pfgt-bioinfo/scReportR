#' scReportR: Single-Cell RNAseq Analysis Report Generator
#'
#' Provides a standardised workflow for single-cell RNAseq analysis.
#' Use [init_sc_report()] to initialise a new project, edit `params.yml`,
#' then render `analysis.qmd` section by section.
#'
#' @keywords internal
"_PACKAGE"

# ── Global imports ────────────────────────────────────────────────────────────
#' @importFrom rlang `%||%`
#' @importFrom dplyr mutate filter select group_by summarise slice_max pull
#'   bind_rows count rename all_of any_of across everything
#' @importFrom tibble tibble rownames_to_column
#' @importFrom ggplot2 ggplot aes geom_violin geom_point geom_col geom_text
#'   geom_hline geom_segment geom_smooth geom_density geom_errorbar
#'   scale_fill_manual scale_color_manual scale_color_viridis_c
#'   scale_y_continuous scale_x_log10 scale_y_log10
#'   facet_wrap facet_grid labs theme theme_bw element_text element_blank
#'   element_rect expansion wrap_plots ggtitle
#' @importFrom patchwork plot_annotation plot_layout
#' @importFrom cli cli_alert_info cli_alert_success cli_alert_warning
#'   cli_alert_danger cli_h1 cli_h2
NULL
