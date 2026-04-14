algorithm_name_to_bnlearn_function <- function(algorithm_name) {
  if (!requireNamespace("bnlearn", quietly = TRUE)) {
    stop("The bnlearn backend requires the bnlearn package.", call. = FALSE)
  }
  switch(
    algorithm_name,
    hillClimbing = bnlearn::hc,
    tabu         = bnlearn::tabu,
    stop("Unsupported bnlearn algorithm name: ", algorithm_name, call. = FALSE)
  )
}

algorithm_name_to_bnlearn_bootstrap_name <- function(algorithm_name) {
  switch(
    algorithm_name,
    hillClimbing = "hc",
    tabu         = "tabu",
    stop("Unsupported bnlearn bootstrap algorithm name: ", algorithm_name, call. = FALSE)
  )
}

fit_bnlearn_networks <- function(model_data, blacklist, configuration) {
  if (!requireNamespace("bnlearn", quietly = TRUE)) {
    stop("The bnlearn backend requires the bnlearn package.", call. = FALSE)
  }
  learned_networks <- list()
  network_scores <- numeric(0)

  for (algorithm_name in configuration$structureLearning$algorithms) {
    learning_function <- algorithm_name_to_bnlearn_function(algorithm_name)
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

learn_bnlearn_network <- function(model_data, blacklist_result, discretizationResult, configuration) {
  if (!requireNamespace("bnlearn", quietly = TRUE)) {
    stop("The bnlearn backend requires the bnlearn package.", call. = FALSE)
  }

  set.seed(configuration$study$randomSeed)
  learned_result <- fit_bnlearn_networks(model_data, blacklist_result$blacklist, configuration)

  bootstrap_algorithm <- algorithm_name_to_bnlearn_bootstrap_name(
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
      modelData          = model_data,
      continuousData     = discretizationResult$continuousAnalysisData,
      discreteData       = discretizationResult$discreteData,
      learnedNetworks    = learned_result$learnedNetworks,
      networkScores      = learned_result$networkScores,
      bootstrapStrength  = bootstrap_strength,
      averagedNetwork    = averaged_network,
      blacklist          = blacklist_result$blacklist,
      variableOrder      = blacklist_result$variableOrder,
      activeExogenousVariables = blacklist_result$activeExogenousVariables,
      measurementGroups  = blacklist_result$measurementGroups,
      measurementColumns = blacklist_result$measurementColumns,
      parameterFit       = parameter_fit,
      threshold          = configuration$structureLearning$averageNetworkThreshold,
      sevtModel          = NULL
    ),
    class = "networkRStructureResult"
  )
}
