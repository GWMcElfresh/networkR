algorithm_name_to_stagedtrees_function <- function(algorithm_name) {
  switch(
    algorithm_name,
    bhc    = stagedtrees::stages_bhc,
    bhcr   = stagedtrees::stages_bhcr,
    bj     = stagedtrees::stages_bj,
    hclust = stagedtrees::stages_hclust,
    kmeans = stagedtrees::stages_kmeans,
    fbhc   = stagedtrees::stages_fbhc,
    stop("Unsupported stagedtrees algorithm: ", algorithm_name, call. = FALSE)
  )
}

fit_stagedtrees_model <- function(model_data, variable_order, configuration) {
  algorithm_name <- configuration$structureLearning$stagedtreesAlgorithm %||% "bhc"
  algorithm_fn <- algorithm_name_to_stagedtrees_function(algorithm_name)
  lambda <- configuration$structureLearning$stagedtreesLambda %||% 0

  full_model <- stagedtrees::full(
    model_data,
    order = variable_order,
    join_unobserved = TRUE,
    lambda = lambda
  )
  primary_model <- algorithm_fn(full_model)
  stagedtrees::sevt_fit(primary_model, data = model_data, lambda = lambda)
}

bnlearn_network_to_matrix <- function(network_object, variable_order) {
  adjacency_matrix <- matrix(
    0L,
    nrow = length(variable_order),
    ncol = length(variable_order),
    dimnames = list(variable_order, variable_order)
  )

  if (is.null(network_object) || !requireNamespace("bnlearn", quietly = TRUE)) {
    return(adjacency_matrix)
  }

  arc_table <- as.data.frame(bnlearn::arcs(network_object), stringsAsFactors = FALSE)
  if (nrow(arc_table) == 0) {
    return(adjacency_matrix)
  }

  for (row_index in seq_len(nrow(arc_table))) {
    from_name <- arc_table$from[[row_index]]
    to_name <- arc_table$to[[row_index]]
    if (from_name %in% variable_order && to_name %in% variable_order) {
      adjacency_matrix[from_name, to_name] <- 1L
    }
  }

  adjacency_matrix
}

# Run one bootstrap replicate: resample rows, fit full model, run algorithm,
# return an adjacency matrix (nodes x nodes, 1 = arc present).
run_bootstrap_replicate <- function(model_data, variable_order, algorithm_fn, lambda) {
  n <- nrow(model_data)
  resampled <- model_data[sample.int(n, n, replace = TRUE), , drop = FALSE]
  # Rebuild factor levels that may have been lost in subset
  resampled <- as.data.frame(lapply(resampled, function(col) {
    if (is.factor(col)) droplevels(col) else col
  }))
  # Restore all original levels so as_adj_matrix gets a consistent node set
  for (v in variable_order) {
    if (is.factor(resampled[[v]]) && is.factor(model_data[[v]])) {
      levels(resampled[[v]]) <- union(levels(resampled[[v]]), levels(model_data[[v]]))
    }
  }
  tryCatch({
    full_model <- stagedtrees::full(resampled, order = variable_order, join_unobserved = TRUE, lambda = lambda)
    fitted_model <- algorithm_fn(full_model)
    stagedtrees::as_adj_matrix(stagedtrees::as_parentslist(fitted_model))
  }, error = function(e) NULL)
}

# Compute arc frequencies over R bootstrap replicates.
# Returns a data frame with columns: from, to, strength (arc frequency in [0,1]).
bootstrap_stagedtrees_network <- function(model_data, variable_order, algorithm_fn, configuration) {
  R <- configuration$structureLearning$bootstrapReplicates
  lambda <- configuration$structureLearning$stagedtreesLambda %||% 0
  n_cores <- configuration$structureLearning$bootstrapCores %||% 1

  replicate_fn <- function(i) {
    run_bootstrap_replicate(model_data, variable_order, algorithm_fn, lambda)
  }

  if (.Platform$OS.type == "windows" || n_cores == 1L) {
    results <- lapply(seq_len(R), replicate_fn)
  } else {
    results <- parallel::mclapply(seq_len(R), replicate_fn, mc.cores = n_cores)
  }

  # Drop failed replicates
  valid_results <- Filter(Negate(is.null), results)
  n_valid <- length(valid_results)
  if (n_valid == 0) {
    stop("All bootstrap replicates failed in learn_stagedtrees_network().", call. = FALSE)
  }

  # Sum adjacency matrices element-wise then divide by n_valid
  freq_matrix <- Reduce(`+`, valid_results) / n_valid

  # Convert accumulated frequencies to from/to/strength long form
  arc_indices <- which(freq_matrix > 0, arr.ind = TRUE)
  if (nrow(arc_indices) == 0) {
    return(data.frame(from = character(0), to = character(0), strength = numeric(0)))
  }

  data.frame(
    from     = rownames(freq_matrix)[arc_indices[, 1]],
    to       = colnames(freq_matrix)[arc_indices[, 2]],
    strength = freq_matrix[arc_indices],
    stringsAsFactors = FALSE
  )
}

# Build the averaged (thresholded) adjacency matrix from bootstrap strengths.
averaged_stagedtrees_network <- function(bootstrap_strength, variable_order, threshold) {
  n <- length(variable_order)
  adj <- matrix(0L, nrow = n, ncol = n, dimnames = list(variable_order, variable_order))
  above_threshold <- bootstrap_strength[bootstrap_strength$strength >= threshold, , drop = FALSE]
  for (i in seq_len(nrow(above_threshold))) {
    from_v <- above_threshold$from[[i]]
    to_v <- above_threshold$to[[i]]
    if (from_v %in% variable_order && to_v %in% variable_order) {
      adj[from_v, to_v] <- 1L
    }
  }
  adj
}

# Main stagedtrees structure learning entry point.
# Returns a networkRStructureResult with an additional $sevtModel slot.
learn_stagedtrees_network <- function(model_data, variable_order, blacklist_result, configuration) {
  threshold <- configuration$structureLearning$averageNetworkThreshold

  set.seed(configuration$study$randomSeed)
  primary_model <- fit_stagedtrees_model(model_data, variable_order, configuration)

  if (requireNamespace("bnlearn", quietly = TRUE)) {
    bnlearn_result <- suppressWarnings(
      learn_bnlearn_network(
        model_data,
        blacklist_result,
        list(continuousAnalysisData = NULL, discreteData = model_data),
        configuration
      )
    )
    bnlearn_result$averagedNetwork <- bnlearn_network_to_matrix(
      bnlearn_result$averagedNetwork,
      variable_order
    )
    bnlearn_result$learnedNetworks$stagedtrees <- primary_model
    bnlearn_result$networkScores[["stagedtrees"]] <- as.numeric(stats::logLik(primary_model))
    bnlearn_result$sevtModel <- primary_model
    return(structure(bnlearn_result, class = "networkRStructureResult"))
  }

  message("bnlearn is not installed; using stagedtrees-native bootstrap estimation for backend='stagedtrees'.")

  algorithm_name <- configuration$structureLearning$stagedtreesAlgorithm %||% "bhc"
  algorithm_fn <- algorithm_name_to_stagedtrees_function(algorithm_name)

  # Bootstrap for arc stability
  set.seed(configuration$study$randomSeed)
  bootstrap_strength <- bootstrap_stagedtrees_network(model_data, variable_order, algorithm_fn, configuration)

  averaged_adj <- averaged_stagedtrees_network(bootstrap_strength, variable_order, threshold)

  structure(
    list(
      modelData          = model_data,
      continuousData     = NULL,  # populated by LearnBayesianNetwork from discretizationResult
      discreteData       = model_data,
      learnedNetworks    = list(primary = primary_model),
      networkScores      = c(primary = logLik(primary_model)),
      bootstrapStrength  = bootstrap_strength,
      averagedNetwork    = averaged_adj,
      blacklist          = blacklist_result$blacklist,
      variableOrder      = variable_order,
      activeExogenousVariables = blacklist_result$activeExogenousVariables,
      measurementGroups  = blacklist_result$measurementGroups,
      measurementColumns = blacklist_result$measurementColumns,
      parameterFit       = NULL,
      threshold          = threshold,
      sevtModel          = primary_model
    ),
    class = "networkRStructureResult"
  )
}
