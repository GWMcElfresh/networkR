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

test_that("ReadPanels infers wide_by_tissue measurement columns when omitted", {
  tb_file <- normalizePath(
    file.path(testthat::test_path("..", ".."), "data", "subjectIdTable_TB.csv"),
    winslash = "/",
    mustWork = TRUE
  )

  configuration <- BuildConfiguration(list(
    data = list(
      dataPath = dirname(tb_file),
      filePattern = "subjectIdTable_TB\\.csv",
      inputFormat = "wide_by_tissue",
      idColumn = "SubjectId",
      tissueColumn = "Tissue",
      rawTimepointColumn = "Timepoint",
      analysisTimepointColumn = "Timepoint",
      groupColumn = "Vaccine",
      protectionColumn = "Probability_Protect",
      measurementColumns = NULL,
      inferMeasurementColumns = TRUE,
      tissueMap = list(
        "Lung-L" = "LungL",
        "Lung-R" = "LungR"
      )
    ),
    variables = list(
      exogenousVariables = c("Vaccine", "Challenge", "Timepoint"),
      stratifierVariables = c("Probability_Protect"),
      measurementGroups = list(
        list(
          groupName = "lunglCells",
          columnPrefixes = c("LungL_"),
          tier = 2
        ),
        list(
          groupName = "lungrCells",
          columnPrefixes = c("LungR_"),
          tier = 2
        )
      )
    )
  ))

  subject_table <- ReadPanels(configuration)

  expect_true(any(grepl("^LungL_TotalMyeloids$", colnames(subject_table))))
  expect_true(any(grepl("^LungR_TotalMyeloids$", colnames(subject_table))))
  expect_false("LungL_cDNA_ID" %in% colnames(subject_table))
  expect_false("LungR_Probability_Protect" %in% colnames(subject_table))
})