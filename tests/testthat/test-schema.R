tb_csv_path <- function() {
  normalizePath(
    file.path(testthat::test_path("..", ".."), "data", "subjectIdTable_TB.csv"),
    winslash = "/",
    mustWork = TRUE
  )
}

test_that("summarize_tabular_schema profiles the TB file and suggests roles", {
  schema_summary <- networkR:::summarize_tabular_schema(tb_csv_path())

  expect_equal(schema_summary$suggestedInputFormat, "wide_by_tissue")
  expect_equal(schema_summary$roleAssignments$subjectId, "SubjectId")
  expect_equal(schema_summary$roleAssignments$sampleId, "cDNA_ID")
  expect_equal(schema_summary$roleAssignments$tissue, "Tissue")
  expect_equal(schema_summary$roleAssignments$group, "Vaccine")
  expect_equal(schema_summary$roleAssignments$outcomeOrProtection, "Probability_Protect")
  expect_true("Challenge" %in% schema_summary$roleAssignments$exogenousNotStratified)
  expect_true("TotalMyeloids" %in% schema_summary$roleAssignments$measurement)
  expect_true("CFU_Homogenate" %in% schema_summary$roleAssignments$unsure)
  expect_true(grepl("columnRoles:", schema_summary$reviewScaffold, fixed = TRUE))
})

test_that("ambiguous numeric metadata columns are routed to unsure with specific rationale", {
  schema_summary <- networkR:::summarize_tabular_schema(tb_csv_path())

  # All three must remain in unsure
  expect_true("CFU_Homogenate" %in% schema_summary$roleAssignments$unsure)
  expect_true("LungPathScore"  %in% schema_summary$roleAssignments$unsure)
  expect_true("Disease_Severity" %in% schema_summary$roleAssignments$unsure)

  # CFU_Homogenate rationale references bacterial burden / CFU terminology
  cfu_issue <- Find(function(issue) identical(issue$column, "CFU_Homogenate"), schema_summary$unresolvedIssues)
  expect_false(is.null(cfu_issue))
  expect_true(grepl("burden|cfu|colony|pathogen", cfu_issue$rationale, ignore.case = TRUE, perl = TRUE))

  # LungPathScore rationale references pathology / score terminology
  path_issue <- Find(function(issue) identical(issue$column, "LungPathScore"), schema_summary$unresolvedIssues)
  expect_false(is.null(path_issue))
  expect_true(grepl("pathology|histology|score|endpoint", path_issue$rationale, ignore.case = TRUE, perl = TRUE))

  # Disease_Severity rationale references severity / clinical / composite
  severity_issue <- Find(function(issue) identical(issue$column, "Disease_Severity"), schema_summary$unresolvedIssues)
  expect_false(is.null(severity_issue))
  expect_true(grepl("severity|clinical|composite|score", severity_issue$rationale, ignore.case = TRUE, perl = TRUE))
})

test_that("parse_role_scaffold round-trips the TB reviewScaffold", {
  schema_summary <- networkR:::summarize_tabular_schema(tb_csv_path())
  scaffold       <- schema_summary$reviewScaffold

  parsed_roles <- networkR:::parse_role_scaffold(scaffold)

  expect_equal(parsed_roles$subjectId, "SubjectId")
  expect_equal(parsed_roles$tissue,    "Tissue")
  expect_equal(parsed_roles$group,     "Vaccine")
  expect_true("TotalMyeloids" %in% parsed_roles$measurement)
  expect_true("CFU_Homogenate" %in% parsed_roles$unsure)
})

test_that("parse_role_scaffold accepts a bare columnRoles block without top-level key", {
  bare_scaffold <- paste0(
    "  subjectId:\n    - \"SubjectId\"\n",
    "  tissue:\n    - \"Tissue\"\n",
    "  measurement:\n    - \"TotalMyeloids\"\n",
    "  unsure:\n    []\n"
  )

  parsed_roles <- networkR:::parse_role_scaffold(bare_scaffold)
  expect_equal(parsed_roles$subjectId, "SubjectId")
  expect_equal(parsed_roles$tissue,    "Tissue")
  expect_true("TotalMyeloids" %in% parsed_roles$measurement)
  expect_equal(length(parsed_roles$unsure), 0L)
})

test_that("resolve_role_assignments returns config-ready pieces", {
  data_frame <- networkR:::read_tabular_data(tb_csv_path())
  schema_summary <- networkR:::summarize_tabular_schema(tb_csv_path())
  resolved_roles <- networkR:::resolve_role_assignments(data_frame, schema_summary$roleAssignments)

  expect_equal(resolved_roles$idColumn, "SubjectId")
  expect_equal(resolved_roles$tissueColumn, "Tissue")
  expect_equal(resolved_roles$protectionColumn, "Probability_Protect")
  expect_true(all(c("Vaccine", "Challenge", "Timepoint") %in% resolved_roles$exogenousVariables))
  expect_true("Probability_Protect" %in% resolved_roles$stratifierVariables)
  expect_true("TotalMyeloids" %in% resolved_roles$measurementColumns)
  expect_false("cDNA_ID" %in% resolved_roles$measurementColumns)
  expect_true("cDNA_ID" %in% resolved_roles$ignoredColumns)
  expect_equal(resolved_roles$tissueMapCandidates[["Lung-L"]], "LungL")
  expect_equal(resolved_roles$tissueMapCandidates[["Lung-R"]], "LungR")
})

test_that("build_minimal_configuration_from_schema generates a compact wide-data config", {
  generated <- networkR:::build_minimal_configuration_from_schema(
    tb_csv_path(),
    study_name = "TB Lung Cell Abundance Network",
    output_directory = "outputs/tb_lung"
  )
  built_configuration <- BuildConfiguration(generated$configuration)

  expect_equal(built_configuration$data$inputFormat, "wide_by_tissue")
  expect_true(isTRUE(built_configuration$data$inferMeasurementColumns))
  expect_equal(built_configuration$data$idColumn, "SubjectId")
  expect_equal(built_configuration$data$tissueColumn, "Tissue")
  expect_equal(built_configuration$data$groupColumn, "Vaccine")
  expect_true("cDNA_ID" %in% built_configuration$data$ignoredColumns)
  expect_true("Probability_Protect" %in% built_configuration$variables$stratifierVariables)
  expect_true("Vaccine" %in% built_configuration$variables$exogenousVariables)
  expect_true(grepl("inferMeasurementColumns:", generated$yamlText, fixed = TRUE))
  expect_true(grepl("ignoredColumns:", generated$yamlText, fixed = TRUE))
  expect_true(grepl("# Overlap guidance:", generated$yamlText, fixed = TRUE))
  expect_true(grepl("Lung-L", generated$yamlText, fixed = TRUE))
})

test_that("BuildConfiguration derives overlap fields when omitted", {
  built_configuration <- BuildConfiguration(list(
    data = list(
      dataPath = "../data",
      idColumn = "SubjectId",
      rawTimepointColumn = "Timepoint",
      analysisTimepointColumn = NULL,
      groupColumn = "Vaccine",
      protectionColumn = "Probability_Protect",
      tissueColumn = "Tissue",
      inputFormat = "wide_by_tissue",
      inferMeasurementColumns = TRUE
    ),
    variables = list(
      exogenousVariables = character(0),
      stratifierVariables = character(0),
      measurementGroups = list(
        list(groupName = "defaultMeasurements", tier = 2)
      )
    )
  ))

  expect_equal(built_configuration$data$analysisTimepointColumn, "Timepoint")
  expect_true(all(c("Vaccine", "Timepoint") %in% built_configuration$variables$exogenousVariables))
  expect_true("Probability_Protect" %in% built_configuration$variables$stratifierVariables)
})