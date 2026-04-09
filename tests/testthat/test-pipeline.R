test_that("RunPipeline produces a pipeline result from CSV panels", {
  panel_directory <- file.path(tempdir(), paste0("networkR_pipeline_", as.integer(Sys.time())))
  dir.create(panel_directory, recursive = TRUE, showWarnings = FALSE)

  liver_kupffer <- c(0.20, 0.25, 0.30, 0.35, 0.28, 0.33, 0.38, 0.43)
  pbmc_monocytes <- c(0.60, 0.56, 0.52, 0.48, 0.58, 0.54, 0.50, 0.46)

  liver_table <- data.frame(
    SubjectId = rep(paste0("S", 1:8), each = 2),
    TimepointLabel2 = "D1",
    Protection = rep(c("Yes", "No", "Yes", "No", "Yes", "No", "Yes", "No"), each = 2),
    GroupName = rep(c("Malaria Preexposed", "Malaria Preexposed", "Malaria Naive", "Malaria Naive", "Malaria Preexposed", "Malaria Preexposed", "Malaria Naive", "Malaria Naive"), each = 2),
    CellType = rep(c("Kupffer Cells", "CD8 T Cells"), times = 8),
    TotalCells = rep(c(20, 80), times = 8),
    TotalCellsForSample = 100,
    Fraction = as.vector(rbind(liver_kupffer, 1 - liver_kupffer))
  )
  pbmc_table <- liver_table
  pbmc_table$CellType <- rep(c("Monocytes", "CD8 T Cells"), times = 8)
  pbmc_table$Fraction <- as.vector(rbind(pbmc_monocytes, 1 - pbmc_monocytes))

  utils::write.csv(liver_table, file.path(panel_directory, "liver_panel.csv"), row.names = FALSE)
  utils::write.csv(pbmc_table, file.path(panel_directory, "pbmc_panel.csv"), row.names = FALSE)

  configuration <- BuildConfiguration(list(
    data = list(
      dataPath = panel_directory,
      filePattern = "[.](csv)$"
    ),
    structureLearning = list(
      algorithms = c("tabu"),
      selectedBootstrapAlgorithm = "tabu",
      bootstrapReplicates = 5,
      averageNetworkThreshold = 0.2
    ),
    imputation = list(
      numberOfImputations = 2,
      maximumIterations = 2
    ),
    stratifications = list(),
    output = list(savePlots = FALSE)
  ))

  pipeline_result <- RunPipeline(configuration)

  expect_s3_class(pipeline_result, "networkRPipelineResult")
  expect_true("fullModel" %in% names(pipeline_result$plotObjects))
})