#' Learn a network structure from discretized measurements
#'
#' Dispatches to the configured backend (`stagedtrees` or `bnlearn`).
#'
#' @param discretizationResult A result from `DiscretizeMeasurements()`.
#' @param configuration A networkR configuration object.
#'
#' @return A `networkRStructureResult` list. When backend is `"stagedtrees"`,
#'   the result also contains a `$sevtModel` slot holding the fitted `sevt`.
LearnBayesianNetwork <- function(discretizationResult, configuration) {
  blacklist_result <- BuildBlacklist(discretizationResult, configuration)
  active_exogenous_variables <- blacklist_result$activeExogenousVariables
  model_columns <- c(active_exogenous_variables, blacklist_result$measurementColumns)
  model_data <- discretizationResult$discreteData[, model_columns, drop = FALSE]

  if (ncol(model_data) < 2) {
    stop("At least two varying model columns are required to learn a network.", call. = FALSE)
  }

  backend <- configuration$structureLearning$backend %||% "stagedtrees"

  if (identical(backend, "stagedtrees")) {
    result <- learn_stagedtrees_network(
      model_data,
      blacklist_result$variableOrder,
      blacklist_result,
      configuration
    )
    result$continuousData <- discretizationResult$continuousAnalysisData
    result
  } else if (identical(backend, "bnlearn")) {
    learn_bnlearn_network(model_data, blacklist_result, discretizationResult, configuration)
  } else {
    stop("Unknown backend: ", backend, call. = FALSE)
  }
}
