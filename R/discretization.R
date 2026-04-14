cut_to_factor <- function(numeric_vector, n_breaks) {
  quantile_breaks <- unique(stats::quantile(numeric_vector, probs = seq(0, 1, length.out = n_breaks + 1), na.rm = TRUE))
  if (length(quantile_breaks) < 2) {
    return(factor(rep("all", length(numeric_vector)), levels = "all"))
  }
  labels <- paste0("q", seq_len(length(quantile_breaks) - 1))
  base::cut(numeric_vector, breaks = quantile_breaks, labels = labels, include.lowest = TRUE)
}

fallback_discretize_vector <- function(numeric_vector, fallback_breaks) {
  if (requireNamespace("bnlearn", quietly = TRUE)) {
    tryCatch(
      bnlearn::discretize(data.frame(value = numeric_vector), method = "quantile", breaks = fallback_breaks)[[1]],
      error = function(error_condition) {
        base::cut(
          numeric_vector,
          breaks = c(-Inf, stats::median(numeric_vector, na.rm = TRUE), Inf),
          labels = c("low", "high"),
          include.lowest = TRUE
        )
      }
    )
  } else {
    cut_to_factor(numeric_vector, fallback_breaks)
  }
}

#' Discretize continuous measurements for Bayesian network learning
#'
#' @param imputationResult An imputation result created by `ImputeMeasurements()`.
#' @param configuration A networkR configuration object.
#'
#' @return A list containing continuous and discrete analysis data.
DiscretizeMeasurements <- function(imputationResult, configuration) {
  if (!inherits(imputationResult, "networkRImputationResult")) {
    stop("imputationResult must be created by ImputeMeasurements().", call. = FALSE)
  }

  continuous_analysis_data <- imputationResult$continuousAnalysisData
  measurement_columns <- imputationResult$measurementColumns
  exogenous_variables <- imputationResult$exogenousVariables
  discrete_data <- continuous_analysis_data

  discrete_data[exogenous_variables] <- lapply(discrete_data[exogenous_variables], as.factor)
  discrete_data[measurement_columns] <- lapply(discrete_data[measurement_columns], function(column_values) {
    if (is.integer(column_values)) as.numeric(column_values) else column_values
  })

  # Resolve method: auto-degrade hartemink to factors for stagedtrees backend
  disc_method <- configuration$discretization$method %||% "factors"
  backend <- configuration$structureLearning$backend %||% "stagedtrees"

  if (disc_method == "hartemink" && backend == "stagedtrees") {
    message("Discretization method 'hartemink' requires bnlearn; switching to 'factors' for stagedtrees backend.")
    disc_method <- "factors"
  }

  fallback_breaks  <- configuration$discretization$fallbackBreaks  %||% 2L
  final_breaks     <- configuration$discretization$finalBreaks      %||% 3L
  initial_breaks   <- configuration$discretization$initialBreaks    %||% 20L

  # ------------------------------------------------------------------ factors
  if (disc_method == "factors") {
    for (col in measurement_columns) {
      discrete_data[[col]] <- cut_to_factor(continuous_analysis_data[[col]], final_breaks)
    }
    fallback_columns  <- character(0)
    hartemink_columns <- character(0)

  # -------------------------------------------------------- hartemink/quantile
  } else {
    if (!requireNamespace("bnlearn", quietly = TRUE)) {
      stop("bnlearn is required for discretization methods 'hartemink' and 'quantile'. ",
           "Install it with install.packages('bnlearn'), or set discretization$method to 'factors'.",
           call. = FALSE)
    }

    can_use_hartemink <- vapply(measurement_columns, function(column_name) {
      tryCatch({
        bnlearn::discretize(
          data.frame(value = continuous_analysis_data[[column_name]]),
          method = "quantile",
          breaks = initial_breaks
        )
        TRUE
      }, error = function(error_condition) FALSE)
    }, logical(1))

    fallback_columns  <- measurement_columns[!can_use_hartemink]
    hartemink_columns <- measurement_columns[can_use_hartemink]

    if (length(fallback_columns) > 0) {
      for (fallback_column in fallback_columns) {
        discrete_data[[fallback_column]] <- fallback_discretize_vector(
          continuous_analysis_data[[fallback_column]],
          fallback_breaks
        )
      }
    }

    if (length(hartemink_columns) >= 2) {
      discrete_data[hartemink_columns] <- bnlearn::discretize(
        continuous_analysis_data[hartemink_columns],
        method = disc_method,
        breaks = final_breaks,
        ibreaks = initial_breaks,
        idisc = "quantile"
      )
    } else if (length(hartemink_columns) == 1) {
      discrete_data[hartemink_columns] <- bnlearn::discretize(
        continuous_analysis_data[hartemink_columns],
        method = "quantile",
        breaks = final_breaks
      )
    }
  }

  non_factor_columns <- measurement_columns[!vapply(discrete_data[measurement_columns], is.factor, logical(1))]
  if (length(non_factor_columns) > 0) {
    stop(
      "Discretization failed for columns: ",
      paste(non_factor_columns, collapse = ", "),
      call. = FALSE
    )
  }

  list(
    subjectTable           = imputationResult$subjectTable,
    imputedData            = imputationResult$imputedData,
    continuousAnalysisData = continuous_analysis_data,
    discreteData           = discrete_data,
    measurementColumns     = measurement_columns,
    exogenousVariables     = exogenous_variables,
    stratifierVariables    = imputationResult$stratifierVariables,
    measurementGroups      = imputationResult$measurementGroups,
    fallbackColumns        = fallback_columns,
    harteminkColumns       = hartemink_columns
  )
}
