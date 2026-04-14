test_that("BuildConfiguration returns a validated configuration", {
  configuration <- BuildConfiguration()

  expect_s3_class(configuration, "networkRConfiguration")
  expect_equal(configuration$study$studyName, "networkR analysis")
  expect_true("GroupName" %in% configuration$variables$exogenousVariables)
})

test_that("ParseConfiguration reads fixture YAML", {
  configuration <- ParseConfiguration(testthat::test_path("fixtures", "minimal_configuration.yml"))

  expect_s3_class(configuration, "networkRConfiguration")
  expect_equal(configuration$study$studyName, "Fixture Study")
  expect_equal(configuration$data$filePattern, "[.](csv)$")
})

test_that("wide_by_tissue configurations can infer measurement columns", {
  configuration <- BuildConfiguration(list(
    data = list(
      dataPath = tempdir(),
      inputFormat = "wide_by_tissue",
      tissueColumn = "Tissue",
      measurementColumns = NULL,
      inferMeasurementColumns = TRUE
    )
  ))

  expect_true(isTRUE(configuration$data$inferMeasurementColumns))
})