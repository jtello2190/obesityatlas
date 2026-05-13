#' Get top marker genes for a Leiden cluster within a tissue
#'
#' Filter a marker table to one tissue and one Leiden cluster, keep statistically
#' significant positive markers, and return the top markers ranked by decreasing
#' average log2 fold change.
#'
#' @param markers A data frame or tibble containing marker results. Must include
#'   the columns `tissue`, `leiden`, `avg_log2FC`, and `p_val_adj`. The table is
#'   expected to have one row per marker feature/gene per tissue-cluster
#'   comparison.
#' @param tissue A single character string giving the tissue to filter, for
#'   example `"liver"` or `"heart"`. The value must be present in
#'   `markers$tissue`.
#' @param cluster A single Leiden cluster identifier. Can be character or
#'   numeric. Internally, both `cluster` and `markers$leiden` are compared as
#'   character values.
#' @param top_n Integer. Maximum number of marker rows to return. Defaults to
#'   `20`.
#' @param min_log2fc Numeric. Minimum `avg_log2FC` required for a marker to be
#'   retained. Defaults to `0.25`.
#'
#' @return A tibble containing up to `top_n` marker rows for the requested tissue
#'   and cluster, filtered to `avg_log2FC >= min_log2fc` and `p_val_adj < 0.05`,
#'   sorted by decreasing `avg_log2FC`.
#'
#' @details
#' This function is intended for marker tables generated from Leiden-cluster
#' differential expression analyses. Positive `avg_log2FC` values are interpreted
#' as enrichment in the target cluster relative to the comparison group. Rows with
#' missing values in `avg_log2FC` or `p_val_adj` are implicitly removed by the
#' filtering conditions.
#'
#' @examples
#' markers <- tibble::tibble(
#'   tissue = c("liver", "liver", "liver", "heart"),
#'   leiden = c("1", "1", "2", "1"),
#'   gene = c("GeneA", "GeneB", "GeneC", "GeneD"),
#'   avg_log2FC = c(1.2, 0.8, 0.4, 1.5),
#'   p_val_adj = c(0.001, 0.02, 0.04, 0.01)
#' )
#'
#' top_markers(markers, tissue = "liver", cluster = "1")
#' top_markers(markers, tissue = "liver", cluster = 1, top_n = 1)
#'
#' @importFrom rlang .data .env
#' @export
top_markers <- function(markers, tissue, cluster,
                        top_n = 20, min_log2fc = 0.25) {
  if (!tissue %in% unique(markers$tissue)) {
    cli::cli_abort("Tissue {.val {tissue}} not found.")
  }
  markers |>
    dplyr::filter(
      .data$tissue == .env$tissue,
      as.character(.data$leiden) == as.character(.env$cluster),
      .data$avg_log2FC >= min_log2fc,
      .data$p_val_adj < 0.05
    ) |>
    dplyr::arrange(dplyr::desc(.data$avg_log2FC)) |>
    dplyr::slice_head(n = top_n)
}

#' Compute cell type composition per tissue and diet
#'
#' @param obs Observation metadata table. Must contain grouping columns in `by`
#'   and a `cell_type` column.
#' @param by Character vector of columns used to group cells before computing
#'   proportions. Defaults to `c("tissue", "diet")`.
#'
#' @returns A tibble with grouping columns, `cell_type`, `n`, and `prop`.
#' @export
#' @importFrom rlang .data
#'
#' @examples
#' obs <- tibble::tibble(
#'   tissue = c("liver", "liver", "heart"),
#'   diet = c("NCD", "HFD", "NCD"),
#'   cell_type = c("A", "A", "B")
#' )
#' cell_composition(obs)
cell_composition <- function(obs, by = c("tissue", "diet")) {
  obs |>
    dplyr::count(dplyr::across(dplyr::all_of(by)), .data$cell_type) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(by))) |>
    dplyr::mutate(prop = .data$n / sum(.data$n)) |>
    dplyr::ungroup()
}

#' Compare HFD vs NCD markers for a cluster
#'
#' @param diet_markers Diet marker tibble. Must contain `tissue`, `leiden`,
#'   `avg_log2FC`, and `p_val_adj`.
#' @param tissue Tissue name to filter.
#' @param cluster Leiden cluster ID to filter.
#' @param top_n Number of top markers to return. Defaults to 20.
#'
#' @returns A tibble of significant diet-associated markers sorted by absolute
#'   log2 fold change.
#' @export
#' @importFrom rlang .data .env
#'
#' @examples
#' diet_markers <- tibble::tibble(
#'   tissue = c("liver", "liver"),
#'   leiden = c("1", "1"),
#'   gene = c("GeneA", "GeneB"),
#'   avg_log2FC = c(1.2, -0.8),
#'   p_val_adj = c(0.01, 0.03)
#' )
#' diet_effect(diet_markers, "liver", "1")
diet_effect <- function(diet_markers, tissue, cluster, top_n = 20) {
  diet_markers |>
    dplyr::filter(
      .data$tissue == .env$tissue,
      .data$leiden == .env$cluster,
      .data$p_val_adj < 0.05
    ) |>
    dplyr::arrange(dplyr::desc(abs(.data$avg_log2FC))) |>
    dplyr::slice_head(n = top_n)
}
