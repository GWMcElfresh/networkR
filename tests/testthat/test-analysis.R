build_fixture_subject_table <- function() {
  data.frame(
    SubjectId = paste0("S", 1:8),
    GroupName = c("Malaria Preexposed", "Malaria Preexposed", "Malaria Naive", "Malaria Naive", "Malaria Preexposed", "Malaria Preexposed", "Malaria Naive", "Malaria Naive"),
    VaccinationStatus = "Vaccinated",
    Timepoint = c("D1", "D1", "D1", "D1", "D2", "D2", "D2", "D2"),
    Protection = c("Yes", "No", "Yes", "No", "Yes", "No", "Yes", "No"),
    Liver_Pct_Kupffer = c(0.20, 0.30, 0.40, 0.45, 0.25, NA, 0.35, 0.50),
    PBMC_Pct_CD8 = c(0.55, 0.60, 0.40, 0.35, 0.52, 0.48, NA, 0.32),
    stringsAsFactors = FALSE
  )
}

test_that("ImputeMeasurements and DiscretizeMeasurements return aligned outputs", {
  configuration <- ParseConfiguration(testthat::test_path("fixtures", "minimal_configuration.yml"))
  subject_table <- build_fixture_subject_table()

  imputation_result <- ImputeMeasurements(subject_table, configuration)
  discretization_result <- DiscretizeMeasurements(imputation_result, configuration)

  expect_true(nrow(imputation_result$imputedData) == nrow(discretization_result$discreteData))
  expect_true(all(vapply(discretization_result$discreteData[discretization_result$measurementColumns], is.factor, logical(1))))
})

test_that("BuildBlacklist adds exogenous and tier constraints", {
  configuration <- ParseConfiguration(testthat::test_path("fixtures", "minimal_configuration.yml"))
  subject_table <- build_fixture_subject_table()
  discretization_result <- DiscretizeMeasurements(ImputeMeasurements(subject_table, configuration), configuration)
  blacklist_result <- BuildBlacklist(discretization_result, configuration)

  expect_true(any(blacklist_result$blacklist$to == "GroupName"))
  expect_true(any(blacklist_result$blacklist$from == "PBMC_Pct_CD8" & blacklist_result$blacklist$to == "Liver_Pct_Kupffer"))
})