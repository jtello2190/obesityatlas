#' Launch the obesity atlas Shiny app
#'
#' Starts the interactive Shiny application bundled with the package.
#' The app explores cell type composition, UMAP embeddings, and marker
#' genes from the mouse multi-tissue obesity atlas.
#'
#' @param data_path Path to the directory containing the atlas CSV files
#'   (`obs.csv`, `umap_atlas.csv`, `leiden_markers.csv`,
#'   `leiden_diet_markers.csv`, `var.csv`, `all_genes.csv`).
#' @param ... Additional arguments passed to [shiny::runApp()].
#'
#' @return Called for its side effect of launching the Shiny app.
#'   Returns invisibly.
#' @export
#'
#' @examples
#' \dontrun{
#'   obesityatlas::run_app("~/data/atlas-csvs")
#' }
run_app <- function(data_path, ...) {
  rlang::check_installed(c("DT", "ggplot2", "shiny"),
                         reason = "to run the Shiny app")

  if (!dir.exists(data_path)) {
    cli::cli_abort("{.arg data_path} does not exist: {.path {data_path}}")
  }

  app_dir <- system.file("app", package = "obesityatlas")
  if (app_dir == "" || !file.exists(file.path(app_dir, "app.R"))) {
    cli::cli_abort(
      "Could not find the app directory. Try reinstalling the package."
    )
  }

  # Pass the data path through an option, since inst/app/app.R
  # runs in its own environment and can't see function arguments
  old_opt <- options(obesityatlas.data_path = normalizePath(data_path))
  on.exit(options(old_opt), add = TRUE)

  shiny::runApp(app_dir, ...)
}
