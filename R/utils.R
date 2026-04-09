assert_has_columns <- function(data_frame, required_columns, object_name = "data_frame") {
  missing_columns <- setdiff(required_columns, colnames(data_frame))
  if (length(missing_columns) > 0) {
    stop(
      object_name,
      " is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
}

utils::globalVariables(c(".data", "measurementValue", "edge_width", "edge_alpha"))

parse_tissue_name <- function(file_path, tissue_map) {
  file_name <- tolower(basename(file_path))
  for (tissue_key in names(tissue_map)) {
    if (grepl(tissue_key, file_name, fixed = TRUE)) {
      return(unname(tissue_map[[tissue_key]]))
    }
  }

  tools::toTitleCase(sub("_.*", "", tools::file_path_sans_ext(basename(file_path))))
}

clean_cell_type <- function(cell_types, cell_type_prefix) {
  cleaned_cell_types <- gsub("[^a-zA-Z0-9]+", "_", cell_types)
  cleaned_cell_types <- gsub("^_|_$", "", cleaned_cell_types)
  paste0(cell_type_prefix, cleaned_cell_types)
}

apply_cell_type_renames <- function(cell_types, tissue_name, rename_rules, cell_type_prefix) {
  if (length(rename_rules) == 0) {
    return(cell_types)
  }

  renamed_cell_types <- cell_types

  for (rename_rule in rename_rules) {
    match_label <- paste0(tissue_name, "_", cell_type_prefix, rename_rule$match)
    replacement_label <- paste0(tissue_name, "_", cell_type_prefix, rename_rule$replace)
    except_tissue <- rename_rule$exceptTissue %||% NULL

    if (!is.null(except_tissue) && identical(tissue_name, except_tissue)) {
      next
    }

    renamed_cell_types[renamed_cell_types == match_label] <- replacement_label
  }

  renamed_cell_types
}

`%||%` <- function(left_value, right_value) {
  if (is.null(left_value)) right_value else left_value
}

find_active_exogenous_variables <- function(data_frame, candidate_names) {
  candidate_names[vapply(candidate_names, function(variable_name) {
    variable_name %in% colnames(data_frame) &&
      nlevels(droplevels(as.factor(data_frame[[variable_name]]))) > 1
  }, logical(1))]
}

find_active_measurement_columns <- function(data_frame, candidate_names) {
  candidate_names[vapply(candidate_names, function(column_name) {
    if (!column_name %in% colnames(data_frame)) {
      return(FALSE)
    }

    column_values <- data_frame[[column_name]]
    if (is.factor(column_values)) {
      nlevels(droplevels(column_values)) > 1
    } else {
      length(unique(stats::na.omit(column_values))) > 1
    }
  }, logical(1))]
}

filter_measurement_groups <- function(measurement_groups, active_measurement_columns) {
  Filter(function(measurement_group) {
    length(measurement_group$columns) > 0
  }, lapply(measurement_groups, function(measurement_group) {
    measurement_group$columns <- intersect(measurement_group$columns, active_measurement_columns)
    measurement_group
  }))
}

match_measurement_groups <- function(column_names, measurement_group_configuration, plotting_configuration) {
  matched_groups <- lapply(measurement_group_configuration, function(measurement_group) {
    matched_columns <- unique(unlist(lapply(measurement_group$columnPrefixes, function(column_prefix) {
      grep(paste0("^", column_prefix), column_names, value = TRUE)
    })))

    list(
      groupName = measurement_group$groupName,
      columns = matched_columns,
      tier = measurement_group$tier,
      fillColour = measurement_group$fillColour %||% plotting_configuration$groupColours$defaultMeasurement$fillColour,
      borderColour = measurement_group$borderColour %||% plotting_configuration$groupColours$defaultMeasurement$borderColour
    )
  })

  matched_columns <- unique(unlist(lapply(matched_groups, `[[`, "columns")))
  unmatched_columns <- setdiff(column_names, matched_columns)

  if (length(unmatched_columns) > 0) {
    matched_groups[[length(matched_groups) + 1]] <- list(
      groupName = "defaultMeasurements",
      columns = unmatched_columns,
      tier = max(c(2, vapply(matched_groups, `[[`, numeric(1), "tier")), na.rm = TRUE),
      fillColour = plotting_configuration$groupColours$defaultMeasurement$fillColour,
      borderColour = plotting_configuration$groupColours$defaultMeasurement$borderColour
    )
  }

  matched_groups
}

prepare_analysis_components <- function(subject_table, configuration) {
  identifier_column <- configuration$data$idColumn
  exogenous_variables <- unique(configuration$variables$exogenousVariables)
  stratifier_variables <- unique(configuration$variables$stratifierVariables %||% character(0))
  excluded_columns <- unique(c(identifier_column, exogenous_variables, stratifier_variables))

  measurement_columns <- colnames(subject_table)[vapply(subject_table, is.numeric, logical(1))]
  measurement_columns <- setdiff(measurement_columns, excluded_columns)

  analysis_columns <- unique(c(exogenous_variables, measurement_columns))
  analysis_data <- subject_table[, intersect(analysis_columns, colnames(subject_table)), drop = FALSE]

  measurement_groups <- match_measurement_groups(
    measurement_columns,
    configuration$variables$measurementGroups,
    configuration$plotting
  )

  list(
    subjectTable = subject_table,
    analysisData = as.data.frame(analysis_data),
    measurementColumns = measurement_columns,
    exogenousVariables = exogenous_variables,
    stratifierVariables = stratifier_variables,
    measurementGroups = measurement_groups
  )
}

apply_filters <- function(data_frame, filters) {
  if (length(filters) == 0) {
    return(seq_len(nrow(data_frame)))
  }

  keep_rows <- rep(TRUE, nrow(data_frame))

  for (filter_definition in filters) {
    column_name <- filter_definition$column
    operator_name <- filter_definition$operator
    filter_value <- filter_definition$value %||% NULL

    if (!column_name %in% colnames(data_frame)) {
      stop("Filter column does not exist: ", column_name, call. = FALSE)
    }

    column_values <- data_frame[[column_name]]

    keep_rows <- keep_rows & switch(
      operator_name,
      equals = column_values == filter_value,
      notEquals = column_values != filter_value,
      "in" = column_values %in% filter_value,
      notMissing = !is.na(column_values),
      missing = is.na(column_values),
      stop("Unsupported filter operator: ", operator_name, call. = FALSE)
    )
  }

  which(keep_rows)
}

make_file_name_fragment <- function(text_value) {
  gsub("[^a-zA-Z0-9]+", "_", tolower(text_value))
}

subset_discretization_result <- function(discretization_result, row_indices) {
  subset_discrete_data <- discretization_result$discreteData[row_indices, , drop = FALSE]
  subset_discrete_data[] <- lapply(subset_discrete_data, function(column_values) {
    if (is.factor(column_values)) droplevels(column_values) else column_values
  })

  list(
    subjectTable = discretization_result$subjectTable[row_indices, , drop = FALSE],
    imputedData = discretization_result$imputedData[row_indices, , drop = FALSE],
    continuousAnalysisData = discretization_result$continuousAnalysisData[row_indices, , drop = FALSE],
    discreteData = subset_discrete_data,
    measurementColumns = discretization_result$measurementColumns,
    exogenousVariables = discretization_result$exogenousVariables,
    stratifierVariables = discretization_result$stratifierVariables,
    measurementGroups = discretization_result$measurementGroups,
    fallbackColumns = discretization_result$fallbackColumns,
    harteminkColumns = discretization_result$harteminkColumns
  )
}

node_group_name <- function(node_name, active_exogenous_variables, measurement_groups) {
  if (node_name %in% active_exogenous_variables) {
    return("exogenous")
  }

  for (measurement_group in measurement_groups) {
    if (node_name %in% measurement_group$columns) {
      return(measurement_group$groupName)
    }
  }

  "unknown"
}

classify_edge_type <- function(node_from, node_to, active_exogenous_variables, measurement_groups) {
  from_group <- node_group_name(node_from, active_exogenous_variables, measurement_groups)
  to_group <- node_group_name(node_to, active_exogenous_variables, measurement_groups)

  if (from_group == "exogenous" && to_group == "exogenous") {
    return("Exogenous to exogenous")
  }

  if (from_group == "exogenous" || to_group == "exogenous") {
    return("Exogenous to measurement")
  }

  if (identical(from_group, to_group)) {
    return(paste("Within", from_group))
  }

  "Cross-compartment"
}