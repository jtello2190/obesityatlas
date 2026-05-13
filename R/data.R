#' Load atlas data from a directory of CSV files
#'
#' @param path Directory containing the atlas CSVs.
#' @return A list with elements `obs`, `markers`, `diet_markers`, `genes`, `var`.
#' @export
#' @importFrom readr read_csv cols col_double
load_atlas <- function(path) {
  if (!dir.exists(path)) {
    cli::cli_abort("Directory {.path {path}} does not exist.")
  }

  required_files <- c(
    "obs.csv",
    "umap_atlas.csv",
    "leiden_markers.csv",
    "leiden_diet_markers.csv",
    "all_genes.csv",
    "var.csv"
  )

  missing_files <- required_files[!file.exists(file.path(path, required_files))]

  if (length(missing_files) > 0) {
    cli::cli_abort(c(
      "Missing required atlas file(s).",
      "x" = "{missing_files}"
    ))
  }

  obs <- readr::read_csv(file.path(path, "obs.csv"),
                         show_col_types = FALSE) |>
    janitor::clean_names()

  umap <- readr::read_csv(
    file.path(path, "umap_atlas.csv"),
    col_names = c("umap_1", "umap_2"),
    col_types = readr::cols(umap_1 = readr::col_double(),
                            umap_2 = readr::col_double()),
    show_col_types = FALSE
  )

  # umap has 128233 rows, obs has 128234 — confirm and bind
  if (nrow(umap) != nrow(obs)) {
    cli::cli_abort(c(
      "UMAP and obs have different numbers of rows.",
      "i" = "obs has {nrow(obs)} rows.",
      "i" = "umap has {nrow(umap)} rows."
    ))
  }

  obs <- dplyr::bind_cols(obs, umap)

  markers <- readr::read_csv(file.path(path, "leiden_markers.csv"),
                             show_col_types = FALSE)
  diet_markers <- readr::read_csv(file.path(path, "leiden_diet_markers.csv"),
                                  show_col_types = FALSE)
  genes <- readr::read_csv(file.path(path, "all_genes.csv"),
                           show_col_types = FALSE)
  var <- readr::read_csv(file.path(path, "var.csv"),
                         show_col_types = FALSE)

  list(
    obs = obs,
    markers = markers,
    diet_markers = diet_markers,
    genes = genes,
    var = var
  )
}
