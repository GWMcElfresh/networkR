#' Impute missing measurement values with MICE
#'
#' @param subjectTable A subject-level table.
#' @param configuration A networkR configuration object.
#'
#' @return A list containing imputed data and analysis metadata.
simple_impute_measurements <- function(analysis_data, measurement_columns) {
  imputed_analysis_data <- analysis_data

  for (measurement_column in measurement_columns) {
    column_values <- imputed_analysis_data[[measurement_column]]
    if (!anyNA(column_values)) {
      next
    }

    replacement_value <- stats::median(column_values, na.rm = TRUE)
    if (is.na(replacement_value)) {
      replacement_value <- 0
    }
    column_values[is.na(column_values)] <- replacement_value
    imputed_analysis_data[[measurement_column]] <- column_values
  }

  imputed_analysis_data
}

#' Impute missing measurement values with MICE
#'
#' @param subjectTable A subject-level table.
#' @param configuration A networkR configuration object.
#'
#' @return A list containing imputed data and analysis metadata.
ImputeMeasurements <- function(subjectTable, configuration) {
  analysis_components <- prepare_analysis_components(subjectTable, configuration)
  analysis_data <- analysis_components$analysisData
  measurement_columns <- analysis_components$measurementColumns

  if (length(measurement_columns) == 0) {
    stop("No numeric measurement columns were detected for imputation.", call. = FALSE)
  }

  continuous_analysis_data <- analysis_data
  imputation_model <- NULL

  if (anyNA(analysis_data[, measurement_columns, drop = FALSE])) {
    imputation_attempt <- tryCatch({
      set.seed(configuration$study$randomSeed)
      initial_imputation <- mice::mice(analysis_data, maxit = 0, printFlag = FALSE)
      predictor_matrix <- initial_imputation$predictorMatrix
      imputation_methods <- initial_imputation$method
      imputation_methods[] <- ""
      imputation_methods[measurement_columns] <- configuration$imputation$method

      excluded_predictors <- intersect(
        configuration$imputation$excludedPredictors,
        colnames(analysis_data)
      )
      predictor_matrix[, excluded_predictors] <- 0

      set.seed(configuration$study$randomSeed)
      imputation_model <- mice::mice(
        analysis_data,
        m = configuration$imputation$numberOfImputations,
        maxit = configuration$imputation$maximumIterations,
        method = imputation_methods,
        predictorMatrix = predictor_matrix,
        printFlag = FALSE
      )

      list(
        continuousAnalysisData = mice::complete(
          imputation_model,
          action = configuration$imputation$completeDataAction
        ),
        imputationModel = imputation_model
      )
    }, error = function(error_condition) {
      warning(
        "MICE imputation failed; using median fallback imputation instead. ",
        conditionMessage(error_condition),
        call. = FALSE
      )
      NULL
    })

    if (is.null(imputation_attempt)) {
      continuous_analysis_data <- simple_impute_measurements(analysis_data, measurement_columns)
    } else {
      continuous_analysis_data <- imputation_attempt$continuousAnalysisData
      imputation_model <- imputation_attempt$imputationModel
    }
  }

  imputed_data <- continuous_analysis_data
  imputed_data[[configuration$data$idColumn]] <- subjectTable[[configuration$data$idColumn]]

  for (stratifier_variable in analysis_components$stratifierVariables) {
    if (stratifier_variable %in% colnames(subjectTable)) {
      imputed_data[[stratifier_variable]] <- subjectTable[[stratifier_variable]]
    }
  }

  remaining_incomplete_rows <- which(!stats::complete.cases(imputed_data[, measurement_columns, drop = FALSE]))
  if (length(remaining_incomplete_rows) > 0) {
    imputed_data <- imputed_data[-remaining_incomplete_rows, , drop = FALSE]
    subjectTable <- subjectTable[-remaining_incomplete_rows, , drop = FALSE]
    continuous_analysis_data <- continuous_analysis_data[-remaining_incomplete_rows, , drop = FALSE]
  }

  structure(
    list(
      subjectTable = subjectTable,
      imputedData = imputed_data,
      continuousAnalysisData = continuous_analysis_data,
      measurementColumns = measurement_columns,
      exogenousVariables = analysis_components$exogenousVariables,
      stratifierVariables = analysis_components$stratifierVariables,
      measurementGroups = analysis_components$measurementGroups,
      imputationModel = imputation_model
    ),
    class = "networkRImputationResult"
  )
}