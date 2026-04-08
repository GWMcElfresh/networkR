test_that("CompareStructures returns delta and edge type columns", {
  first_structure <- list(
    bootstrapStrength = data.frame(from = "PBMC_Pct_CD8", to = "Liver_Pct_Kupffer", strength = 0.6, stringsAsFactors = FALSE),
    continuousData = data.frame(PBMC_Pct_CD8 = c(1, 2, 3, 4, 5), Liver_Pct_Kupffer = c(2, 4, 6, 8, 10)),
    activeExogenousVariables = c("GroupName"),
    measurementGroups = list(
      list(groupName = "pbmcCells", columns = "PBMC_Pct_CD8", tier = 2, fillColour = "lightblue", borderColour = "steelblue"),
      list(groupName = "liverCells", columns = "Liver_Pct_Kupffer", tier = 3, fillColour = "lightgreen", borderColour = "darkgreen")
    )
  )
  second_structure <- list(
    bootstrapStrength = data.frame(from = "PBMC_Pct_CD8", to = "Liver_Pct_Kupffer", strength = 0.2, stringsAsFactors = FALSE),
    continuousData = data.frame(PBMC_Pct_CD8 = c(1, 2, 3, 4, 5), Liver_Pct_Kupffer = c(5, 4, 3, 2, 1)),
    activeExogenousVariables = c("GroupName"),
    measurementGroups = first_structure$measurementGroups
  )

  comparison_result <- CompareStructures(first_structure, second_structure, "First", "Second")

  expect_true(all(c("delta", "absoluteDelta", "edgeType") %in% colnames(comparison_result$comparisonTable)))
  expect_equal(nrow(comparison_result$comparisonTable), 1)
})