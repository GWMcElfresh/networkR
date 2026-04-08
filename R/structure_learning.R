algorithm_name_to_function <- function(algorithm_name) {
  switch(
    algorithm_name,
    hillClimbing = bnlearn::hc,
    tabu = bnlearn::tabu,
    stop("Unsupported algorithm name: ", algorithm_name, call. = FALSE)
  )
}

algorithm_name_to_bootstrap_name <- function(algorithm_name) {
  switch(
    algorithm_name,
    hillClimbing = "hc",
    tabu = "tabu",
    stop("Unsupported bootstrap algorithm name: ", algorithm_name, call. = FALSE)
  )
}

fit_requested_networks <- function(model_data, blacklist, configuration) {
  learned_networks <- list()
  network_scores <- numeric(0)

  for (algorithm_name in configuration$structureLearning$algorithms) {
    learning_function <- algorithm_name_to_function(algorithm_name)
    learned_network <- learning_function(
      model_data,
      blacklist = blacklist,
      score = configuration$structureLearning$score
    )
    learned_networks[[algorithm_name]] <- learned_network
    network_scores[[algorithm_name]] <- bnlearn::score(
      learned_network,
      data = model_data,
      type = configuration$structureLearning$score
    )
  }

  list(learnedNetworks = learned_networks, networkScores = network_scores)
}

#' Learn a Bayesian network from discretized measurements
#'
#' @param discretizationResult A result from `DiscretizeMeasurements()`.
#' @param configuration A networkR configuration object.
#'
#' @return A list containing learned networks, bootstrap strength, and metadata.
LearnBayesianNetwork <- function(discretizationResult, configuration) {
  blacklist_result <- BuildBlacklist(discretizationResult, configuration)
  active_exogenous_variables <- blacklist_result$activeExogenousVariables
  model_columns <- c(active_exogenous_variables, discretizationResult$measurementColumns)
  model_data <- discretizationResult$discreteData[, model_columns, drop = FALSE]

  if (ncol(model_data) < 2) {
    stop("At least two varying model columns are required to learn a network.", call. = FALSE)
  }

  set.seed(configuration$study$randomSeed)
  learned_result <- fit_requested_networks(model_data, blacklist_result$blacklist, configuration)

  bootstrap_algorithm <- algorithm_name_to_bootstrap_name(
    configuration$structureLearning$selectedBootstrapAlgorithm
  )

  set.seed(configuration$study$randomSeed)
  bootstrap_strength <- bnlearn::boot.strength(
    data = model_data,
    R = configuration$structureLearning$bootstrapReplicates,
    algorithm = bootstrap_algorithm,
    algorithm.args = list(
      blacklist = blacklist_result$blacklist,
      score = configuration$structureLearning$score
    )
  )

  averaged_network <- bnlearn::averaged.network(
    bootstrap_strength,
    threshold = configuration$structureLearning$averageNetworkThreshold
  )

  primary_algorithm <- configuration$structureLearning$algorithms[[1]]
  parameter_fit <- NULL
  if (length(learned_result$learnedNetworks) > 0) {
    parameter_fit <- bnlearn::bn.fit(
      learned_result$learnedNetworks[[primary_algorithm]],
      data = model_data,
      method = configuration$structureLearning$parameterMethod,
      iss = configuration$structureLearning$imaginarySampleSize
    )
  }

  structure(
    list(
      modelData = model_data,
      continuousData = discretizationResult$continuousAnalysisData,
      discreteData = discretizationResult$discreteData,
      learnedNetworks = learned_result$learnedNetworks,
      networkScores = learned_result$networkScores,
      bootstrapStrength = bootstrap_strength,
      averagedNetwork = averaged_network,
      blacklist = blacklist_result$blacklist,
      activeExogenousVariables = active_exogenous_variables,
      measurementGroups = discretizationResult$measurementGroups,
      measurementColumns = discretizationResult$measurementColumns,
      parameterFit = parameter_fit,
      threshold = configuration$structureLearning$averageNetworkThreshold
    ),
    class = "networkRStructureResult"
  )
}