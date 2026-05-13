test_that("top_markers filters and sorts correctly", {
  markers <- tibble::tibble(
    tissue = c("Adipose", "Adipose", "Adipose", "Liver"),
    leiden = c("0", "0", "0", "0"),
    geneID = c("Gene1", "Gene2", "Gene3", "Gene4"),
    avg_log2FC = c(2.5, 0.1, 1.0, 3.0),
    p_val_adj = c(0.01, 0.01, 0.01, 0.01),
    pct.1 = c(0.9, 0.5, 0.7, 0.9),
    pct.2 = c(0.1, 0.4, 0.3, 0.1),
    highly_variable = c(TRUE, TRUE, TRUE, TRUE)
  )

  result <- top_markers(markers, "Adipose", "0", top_n = 5, min_log2fc = 0.5)

  expect_equal(nrow(result), 2)
  expect_equal(result$geneID, c("Gene1", "Gene3"))
})

test_that("top_markers errors on unknown tissue", {
  markers <- tibble::tibble(tissue = "Adipose", leiden = "0",
                            avg_log2FC = 1, p_val_adj = 0.01)
  expect_error(top_markers(markers, "Brain", "0"))
})
