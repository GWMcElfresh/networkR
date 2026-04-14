quantile_factor_fixture <- function(values, n_breaks = 3) {
  breaks <- unique(stats::quantile(values, probs = seq(0, 1, length.out = n_breaks + 1), na.rm = TRUE))
  labels <- paste0("q", seq_len(length(breaks) - 1))
  factor(base::cut(values, breaks = breaks, labels = labels, include.lowest = TRUE), levels = labels)
}

build_backend_fixture_configuration <- function(backend, bootstrap_replicates = 60) {
  BuildConfiguration(list(
    study = list(
      studyName = "Backend parity fixture",
      randomSeed = 17
    ),
    data = list(
      dataPath = ".",
      rawTimepointColumn = "Timepoint",
      analysisTimepointColumn = "Timepoint",
      groupColumn = "GroupName"
    ),
    variables = list(
      exogenousVariables = c("GroupName", "Timepoint"),
      stratifierVariables = character(0),
      measurementGroups = list(
        list(groupName = "tier2", columnPrefixes = c("Tier2_"), tier = 2, fillColour = "grey85", borderColour = "grey25"),
        list(groupName = "tier3", columnPrefixes = c("Tier3_"), tier = 3, fillColour = "grey85", borderColour = "grey25"),
        list(groupName = "tier4", columnPrefixes = c("Tier4_"), tier = 4, fillColour = "grey85", borderColour = "grey25"),
        list(groupName = "tier5", columnPrefixes = c("Tier5_"), tier = 5, fillColour = "grey85", borderColour = "grey25")
      )
    ),
    structureLearning = list(
      backend = backend,
      algorithms = c("tabu"),
      selectedBootstrapAlgorithm = "tabu",
      stagedtreesAlgorithm = "bhc",
      bootstrapReplicates = bootstrap_replicates,
      bootstrapCores = 1,
      averageNetworkThreshold = 0.2
    ),
    output = list(savePlots = FALSE)
  ))
}

build_backend_parity_fixture <- function() {
  set.seed(123)

  n_rows <- 240
  group_name <- factor(rep(c("A", "B"), each = n_rows / 2))
  timepoint <- factor(rep(rep(c("T1", "T2"), each = n_rows / 4), 2))

  tier2_signal <- ifelse(group_name == "A", -1.5, 1.5) + stats::rnorm(n_rows, sd = 0.12)
  tier3_signal <- tier2_signal + stats::rnorm(n_rows, sd = 0.10)
  tier4_signal <- tier3_signal + ifelse(timepoint == "T1", -0.8, 0.8) + stats::rnorm(n_rows, sd = 0.10)
  tier5_signal <- tier4_signal + stats::rnorm(n_rows, sd = 0.08)

  continuous_data <- data.frame(
    GroupName = group_name,
    Timepoint = timepoint,
    Tier2_signal = tier2_signal,
    Tier3_signal = tier3_signal,
    Tier4_signal = tier4_signal,
    Tier5_signal = tier5_signal,
    stringsAsFactors = FALSE
  )

  discrete_data <- data.frame(
    GroupName = factor(group_name),
    Timepoint = factor(timepoint),
    Tier2_signal = quantile_factor_fixture(tier2_signal),
    Tier3_signal = quantile_factor_fixture(tier3_signal),
    Tier4_signal = quantile_factor_fixture(tier4_signal),
    Tier5_signal = quantile_factor_fixture(tier5_signal),
    stringsAsFactors = FALSE
  )

  measurement_groups <- list(
    list(groupName = "tier2", columns = "Tier2_signal", tier = 2, fillColour = "grey85", borderColour = "grey25"),
    list(groupName = "tier3", columns = "Tier3_signal", tier = 3, fillColour = "grey85", borderColour = "grey25"),
    list(groupName = "tier4", columns = "Tier4_signal", tier = 4, fillColour = "grey85", borderColour = "grey25"),
    list(groupName = "tier5", columns = "Tier5_signal", tier = 5, fillColour = "grey85", borderColour = "grey25")
  )

  list(
    discretizationResult = list(
      subjectTable = continuous_data,
      imputedData = continuous_data,
      continuousAnalysisData = continuous_data,
      discreteData = discrete_data,
      measurementColumns = c("Tier2_signal", "Tier3_signal", "Tier4_signal", "Tier5_signal"),
      exogenousVariables = c("GroupName", "Timepoint"),
      stratifierVariables = character(0),
      measurementGroups = measurement_groups,
      fallbackColumns = character(0),
      harteminkColumns = character(0)
    ),
    stagedtreesConfiguration = build_backend_fixture_configuration("stagedtrees"),
    bnlearnConfiguration = build_backend_fixture_configuration("bnlearn")
  )
}

test_that("CompareBackendBootstrapStrengths zero-fills missing directed edges", {
  first_structure <- list(
    bootstrapStrength = data.frame(
      from = c("A", "B"),
      to = c("B", "C"),
      strength = c(0.50, 0.40),
      stringsAsFactors = FALSE
    )
  )
  second_structure <- list(
    bootstrapStrength = data.frame(
      from = c("A"),
      to = c("B"),
      strength = c(0.54),
      stringsAsFactors = FALSE
    )
  )

  parity_result <- CompareBackendBootstrapStrengths(first_structure, second_structure, tolerance = 0.05)

  expect_equal(parity_result$summary$edgeCount, 2)
  expect_equal(parity_result$summary$exceedanceCount, 1)
  expect_equal(
    parity_result$edgeTable$secondStrength[parity_result$edgeTable$edge == "B -> C"],
    0
  )
})

test_that("stagedtrees and bnlearn bootstrap strengths stay within parity tolerance on constrained synthetic data", {
  skip_if_not_installed("bnlearn")

  fixture <- build_backend_parity_fixture()

  stagedtrees_result <- LearnBayesianNetwork(
    fixture$discretizationResult,
    fixture$stagedtreesConfiguration
  )
  bnlearn_result <- LearnBayesianNetwork(
    fixture$discretizationResult,
    fixture$bnlearnConfiguration
  )

  parity_result <- CompareBackendBootstrapStrengths(
    stagedtrees_result,
    bnlearn_result,
    firstLabel = "stagedtrees",
    secondLabel = "bnlearn",
    tolerance = 0.05
  )

  expect_gt(parity_result$summary$edgeCount, 0)
  expect_lte(parity_result$summary$exceedanceRate, 0.05)
})

test_that("stagedtrees plotting helpers return ggplot objects", {
  fixture <- build_backend_parity_fixture()

  stagedtrees_result <- LearnBayesianNetwork(
    fixture$discretizationResult,
    fixture$stagedtreesConfiguration
  )

  event_tree_plot <- PlotEventTree(stagedtrees_result, fixture$stagedtreesConfiguration)
  stage_heatmap_plot <- PlotStageSplitHeatmap(stagedtrees_result, fixture$stagedtreesConfiguration)

  expect_s3_class(event_tree_plot, "ggplot")
  expect_s3_class(stage_heatmap_plot, "ggplot")

  layer_geoms <- vapply(event_tree_plot$layers, function(layer) {
    class(layer$geom)[[1]]
  }, character(1))
  expect_true("GeomLabelRepel" %in% layer_geoms)

  y_scale <- event_tree_plot$scales$get_scales("y")
  expect_false(is.null(y_scale))
  expect_true(is.function(y_scale$labels) || length(y_scale$labels) > 0)
})