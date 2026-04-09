create_panel_file <- function(file_path, tissue_label) {
  panel_table <- data.frame(
    SubjectId = rep(paste0("S", 1:8), each = 2),
    TimepointLabel2 = "D1",
    Protection = rep(c("Yes", "No", "Yes", "No", "Yes", "No", "Yes", "No"), each = 2),
    GroupName = rep(c("Malaria Preexposed", "Malaria Preexposed", "Malaria Naive", "Malaria Naive", "Malaria Preexposed", "Malaria Preexposed", "Malaria Naive", "Malaria Naive"), each = 2),
    CellType = rep(c("Kupffer Cells", "CD8 T Cells"), times = 8),
    TotalCells = rep(c(25, 75), times = 8),
    TotalCellsForSample = 100,
    Fraction = rep(c(0.25, 0.75), times = 8)
  )

  utils::write.csv(panel_table, file_path, row.names = FALSE)
}

test_that("ReadPanels reads and joins panel tables", {
  panel_directory <- file.path(tempdir(), paste0("networkR_test_", as.integer(Sys.time())))
  dir.create(panel_directory, recursive = TRUE, showWarnings = FALSE)
  create_panel_file(file.path(panel_directory, "liver_panel.csv"), "Liver")
  create_panel_file(file.path(panel_directory, "pbmc_panel.csv"), "PBMC")

  configuration <- BuildConfiguration(list(
    data = list(
      dataPath = panel_directory,
      filePattern = "[.](csv)$"
    )
  ))

  subject_table <- ReadPanels(configuration)
  subject_table <- AddDerivedVariables(subject_table, configuration)

  expect_equal(nrow(subject_table), 8)
  expect_true("VaccinationStatus" %in% colnames(subject_table))
  expect_false("Tissue" %in% colnames(subject_table))
  expect_true(any(grepl("^Liver_", colnames(subject_table))))
  expect_true(any(grepl("^PBMC_", colnames(subject_table))))
  expect_false(anyNA(subject_table[, grep("^(Liver|PBMC)_", colnames(subject_table), value = TRUE), drop = FALSE]))
})