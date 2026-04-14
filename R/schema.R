normalize_schema_name <- function(column_name) {
  gsub("[^a-z0-9]+", "", tolower(column_name))
}

is_numericish <- function(column_values) {
  is.numeric(column_values) || is.integer(column_values)
}

first_or_null <- function(values) {
  if (length(values) == 0) {
    return(NULL)
  }

  values[[1]]
}

schema_safe_scalar <- function(value) {
  if (length(value) == 0L || all(is.na(value)) || any(is.infinite(value))) {
    return(NA)
  }

  unname(value[[1]])
}

sanitize_tabular_columns <- function(data_frame) {
  valid_name_mask <- nchar(trimws(colnames(data_frame))) > 0 &
    !grepl("^\\.\\.\\.\\d+$", colnames(data_frame))

  data_frame[, valid_name_mask, drop = FALSE]
}

read_tabular_data <- function(file_path) {
  file_extension <- tolower(tools::file_ext(file_path))

  if (identical(file_extension, "csv")) {
    data_frame <- utils::read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (file_extension %in% c("xls", "xlsx")) {
    data_frame <- as.data.frame(readxl::read_excel(file_path), check.names = FALSE)
  } else {
    stop("Unsupported tabular file extension: ", file_extension, call. = FALSE)
  }

  sanitize_tabular_columns(data_frame)
}

suggest_input_format <- function(data_frame) {
  normalized_names <- vapply(colnames(data_frame), normalize_schema_name, character(1))
  has_cell_type <- any(normalized_names %in% c("celltype", "celltypes"))
  has_fraction <- any(normalized_names %in% c("fraction", "fractions"))
  has_tissue <- any(normalized_names %in% c("tissue", "tissues"))

  numeric_columns <- vapply(data_frame, is_numericish, logical(1))

  if (has_cell_type && has_fraction) {
    return("panel")
  }

  if (has_tissue && sum(numeric_columns) > 5) {
    return("wide_by_tissue")
  }

  "panel"
}

summarize_schema_column <- function(name, column, max_unique_values) {
  non_missing_values <- column[!is.na(column)]
  unique_values <- unique(as.character(non_missing_values))

  column_summary <- list(
    name = name,
    class = paste(class(column), collapse = "/"),
    nMissing = sum(is.na(column)),
    nUnique = length(unique_values),
    exampleValues = utils::head(unique_values, max_unique_values)
  )

  if (is_numericish(column)) {
    column_summary$numericSummary <- list(
      min = schema_safe_scalar(min(column, na.rm = TRUE)),
      median = schema_safe_scalar(stats::median(column, na.rm = TRUE)),
      max = schema_safe_scalar(max(column, na.rm = TRUE))
    )
  }

  column_summary
}

empty_role_assignments <- function() {
  stats::setNames(
    replicate(12, character(0), simplify = FALSE),
    c(
      "subjectId",
      "sampleId",
      "tissue",
      "rawTimepoint",
      "analysisTimepoint",
      "group",
      "outcomeOrProtection",
      "exogenousNotStratified",
      "stratifiedBy",
      "measurement",
      "ignored",
      "unsure"
    )
  )
}

role_assignment <- function(role, confidence, rationale) {
  list(
    role = role,
    confidence = confidence,
    rationale = rationale
  )
}

select_named_column <- function(column_summaries, patterns) {
  normalized_names <- vapply(column_summaries, function(column_summary) {
    normalize_schema_name(column_summary$name)
  }, character(1))

  for (pattern in patterns) {
    matched_index <- which(grepl(pattern, normalized_names, perl = TRUE))
    if (length(matched_index) > 0) {
      return(column_summaries[[matched_index[[1]]]]$name)
    }
  }

  NULL
}

suggest_column_roles <- function(data_frame, max_unique_values = 10L) {
  column_summaries <- lapply(colnames(data_frame), function(column_name) {
    summarize_schema_column(column_name, data_frame[[column_name]], max_unique_values)
  })

  role_assignments <- empty_role_assignments()
  unresolved_issues <- list()
  column_lookup <- stats::setNames(vector("list", length(column_summaries)), colnames(data_frame))

  subject_column <- select_named_column(column_summaries, c("^subjectid$", "^subject$", "subject.*id", "animal.*id"))
  sample_column <- select_named_column(column_summaries, c("^cdnaid$", "sample.*id", "specimen.*id", "aliquot.*id"))
  tissue_column <- select_named_column(column_summaries, c("^tissue$"))
  timepoint_column <- select_named_column(column_summaries, c("^timepoint$", "time.*point", "^day$", "collectionday"))
  group_column <- select_named_column(column_summaries, c("^group$", "groupname", "^vaccine$", "treatment", "^arm$", "cohort"))
  outcome_column <- select_named_column(column_summaries, c("protect", "outcome", "response"))

  if (!is.null(subject_column)) {
    role_assignments$subjectId <- subject_column
  }
  if (!is.null(sample_column)) {
    role_assignments$sampleId <- sample_column
  }
  if (!is.null(tissue_column)) {
    role_assignments$tissue <- tissue_column
  }
  if (!is.null(timepoint_column)) {
    role_assignments$rawTimepoint <- timepoint_column
    role_assignments$analysisTimepoint <- timepoint_column
  }
  if (!is.null(group_column)) {
    role_assignments$group <- group_column
  }
  if (!is.null(outcome_column)) {
    role_assignments$outcomeOrProtection <- outcome_column
  }

  reserved_columns <- unique(unlist(role_assignments, use.names = FALSE))

  for (column_summary in column_summaries) {
    column_name <- column_summary$name
    normalized_name <- normalize_schema_name(column_name)
    column_values <- data_frame[[column_name]]

    assignment <- NULL

    if (identical(column_name, subject_column)) {
      assignment <- role_assignment("subjectId", "high", "Column name matches a subject identifier pattern.")
    } else if (identical(column_name, sample_column)) {
      assignment <- role_assignment("sampleId", "high", "Column name suggests a per-measurement sample identifier.")
    } else if (identical(column_name, tissue_column)) {
      assignment <- role_assignment("tissue", "high", "Column name matches the tissue axis for a wide-by-tissue table.")
    } else if (identical(column_name, timepoint_column)) {
      assignment <- role_assignment("rawTimepoint", "high", "Column name matches a timepoint axis and can also serve as analysisTimepoint.")
    } else if (identical(column_name, group_column)) {
      assignment <- role_assignment("group", "high", "Column name matches a treatment or study-group variable.")
    } else if (identical(column_name, outcome_column)) {
      assignment <- role_assignment("outcomeOrProtection", "high", "Column name suggests an outcome or protection variable.")
    } else if (!is_numericish(column_values) || is.logical(column_values)) {
      if (column_summary$nUnique <= 1) {
        assignment <- role_assignment("ignored", "medium", "Column has only one observed level and will not influence the analysis.")
      } else if (column_summary$nUnique >= nrow(data_frame) * 0.8) {
        assignment <- role_assignment("unsure", "low", "Column behaves like an identifier but its name does not match a known identifier pattern.")
      } else {
        assignment <- role_assignment("exogenousNotStratified", "medium", "Categorical column is a plausible design variable rather than a measurement.")
      }
    } else {
      if (column_summary$nUnique <= 1) {
        assignment <- role_assignment("ignored", "medium", "Numeric column is constant and can be ignored.")
      } else if (grepl("protect|outcome|response", normalized_name, perl = TRUE)) {
        assignment <- role_assignment("outcomeOrProtection", "medium", "Numeric column name suggests an outcome or protection signal.")
        role_assignments$outcomeOrProtection <- unique(c(role_assignments$outcomeOrProtection, column_name))
      } else if (grepl("^cfu|cfu[0-9_]|colonyforming|cfucount|cfuhom|pathogenload|bacterialload|baccount", normalized_name, perl = TRUE)) {
        assignment <- role_assignment("unsure", "medium", paste(
          "Numeric column looks like a bacterial burden or pathogen load endpoint (e.g., colony-forming unit counts).",
          "Move to 'stratifiedBy' to split analyses by load, 'outcomeOrProtection' if this is the primary endpoint,",
          "or 'ignored' to exclude it from the network."
        ))
      } else if (grepl("severity", normalized_name, perl = TRUE)) {
        assignment <- role_assignment("unsure", "medium", paste(
          "Numeric column looks like a clinical severity or composite disease scoring variable.",
          "Move to 'stratifiedBy', 'outcomeOrProtection', 'measurement', or 'ignored' as appropriate."
        ))
      } else if (grepl("score|burden|path", normalized_name, perl = TRUE)) {
        assignment <- role_assignment("unsure", "medium", paste(
          "Numeric column looks like a pathology score, histology grade, or analysis endpoint",
          "rather than a flow cytometry cell measurement.",
          "Move to 'stratifiedBy', 'measurement', or 'ignored' as appropriate."
        ))
      } else if (column_summary$nUnique <= 5) {
        assignment <- role_assignment("unsure", "low", "Numeric column has very low cardinality and may be better treated as a stratifier or metadata.")
      } else {
        assignment <- role_assignment("measurement", "high", "Numeric column has enough variation to be treated as a candidate measurement.")
      }
    }

    if (!(column_name %in% reserved_columns) && assignment$role %in% names(role_assignments)) {
      role_assignments[[assignment$role]] <- unique(c(role_assignments[[assignment$role]], column_name))
    }

    if (identical(assignment$role, "unsure")) {
      unresolved_issues[[length(unresolved_issues) + 1]] <- list(
        column = column_name,
        rationale = assignment$rationale
      )
    }

    column_summary$roleSuggestion <- assignment
    column_lookup[[column_name]] <- column_summary
  }

  required_roles <- c("subjectId", "tissue", "rawTimepoint", "group")
  for (required_role in required_roles) {
    if (length(role_assignments[[required_role]]) == 0) {
      unresolved_issues[[length(unresolved_issues) + 1]] <- list(
        column = required_role,
        rationale = paste("No column was confidently assigned to role:", required_role)
      )
    }
  }

  list(
    columns = unname(column_lookup),
    roleAssignments = role_assignments,
    unresolvedIssues = unresolved_issues
  )
}

infer_wide_measurement_columns <- function(data_frame, configuration = NULL) {
  explicit_measurements <- configuration$data$measurementColumns %||% NULL
  if (!is.null(explicit_measurements) && length(explicit_measurements) > 0) {
    return(explicit_measurements)
  }

  role_suggestions <- suggest_column_roles(data_frame)
  inferred_measurements <- role_suggestions$roleAssignments$measurement

  if (length(inferred_measurements) > 0) {
    return(inferred_measurements)
  }

  excluded_columns <- unique(c(
    configuration$data$idColumn,
    configuration$data$tissueColumn,
    configuration$data$rawTimepointColumn,
    configuration$data$analysisTimepointColumn,
    configuration$data$groupColumn,
    configuration$data$protectionColumn,
    configuration$variables$exogenousVariables %||% character(0),
    configuration$variables$stratifierVariables %||% character(0)
  ))

  numeric_columns <- colnames(data_frame)[vapply(data_frame, is_numericish, logical(1))]
  setdiff(numeric_columns, excluded_columns)
}

resolve_role_assignments <- function(data_frame, role_assignments) {
  tissue_column <- first_or_null(role_assignments$tissue)
  tissue_map_candidates <- list()

  if (!is.null(tissue_column) && tissue_column %in% colnames(data_frame)) {
    tissue_values <- unique(as.character(stats::na.omit(data_frame[[tissue_column]])))
    tissue_map_candidates <- stats::setNames(
      as.list(gsub("[^A-Za-z0-9]", "", tissue_values)),
      tissue_values
    )
  }

  list(
    idColumn = first_or_null(role_assignments$subjectId),
    sampleIdColumns = role_assignments$sampleId,
    tissueColumn = tissue_column,
    rawTimepointColumn = first_or_null(role_assignments$rawTimepoint),
    analysisTimepointColumn = first_or_null(role_assignments$analysisTimepoint),
    groupColumn = first_or_null(role_assignments$group),
    protectionColumn = first_or_null(role_assignments$outcomeOrProtection),
    exogenousVariables = unique(c(
      first_or_null(role_assignments$group),
      first_or_null(role_assignments$rawTimepoint),
      role_assignments$exogenousNotStratified
    )),
    stratifierVariables = unique(c(
      role_assignments$stratifiedBy,
      first_or_null(role_assignments$outcomeOrProtection)
    )),
    measurementColumns = role_assignments$measurement,
    ignoredColumns = unique(c(role_assignments$ignored, role_assignments$sampleId)),
    tissueMapCandidates = tissue_map_candidates,
    unresolvedColumns = role_assignments$unsure
  )
}

build_column_role_scaffold <- function(schema_summary) {
  role_assignments <- schema_summary$roleAssignments
  column_rationales <- stats::setNames(
    vapply(schema_summary$columns, function(column_summary) {
      column_summary$roleSuggestion$rationale
    }, character(1)),
    vapply(schema_summary$columns, `[[`, character(1), "name")
  )

  scaffold_lines <- c(
    "# Review these column-role assignments.",
    "# Move column names between sections before generating YAML.",
    "columnRoles:"
  )

  for (role_name in names(role_assignments)) {
    scaffold_lines <- c(scaffold_lines, paste0("  ", role_name, ":"))
    role_columns <- role_assignments[[role_name]]

    if (length(role_columns) == 0) {
      scaffold_lines <- c(scaffold_lines, "    []")
      next
    }

    for (column_name in role_columns) {
      rationale <- column_rationales[[column_name]] %||% ""
      if (nzchar(rationale)) {
        scaffold_lines <- c(
          scaffold_lines,
          paste0("    - \"", column_name, "\" # ", rationale)
        )
      } else {
        scaffold_lines <- c(scaffold_lines, paste0("    - \"", column_name, "\""))
      }
    }
  }

  paste(scaffold_lines, collapse = "\n")
}

summarize_tabular_schema <- function(file_path, max_unique_values = 10L) {
  data_frame <- read_tabular_data(file_path)
  suggested_format <- suggest_input_format(data_frame)
  role_suggestions <- suggest_column_roles(data_frame, max_unique_values = max_unique_values)
  resolved_roles <- resolve_role_assignments(data_frame, role_suggestions$roleAssignments)

  list(
    filePath = file_path,
    nRows = nrow(data_frame),
    nCols = ncol(data_frame),
    suggestedInputFormat = suggested_format,
    columns = role_suggestions$columns,
    roleAssignments = role_suggestions$roleAssignments,
    unresolvedIssues = role_suggestions$unresolvedIssues,
    resolvedRoles = resolved_roles,
    reviewScaffold = build_column_role_scaffold(list(
      columns = role_suggestions$columns,
      roleAssignments = role_suggestions$roleAssignments
    ))
  )
}

parse_role_scaffold <- function(scaffold_yaml) {
  parsed <- yaml::yaml.load(scaffold_yaml)

  roles <- if (!is.null(parsed$columnRoles)) parsed$columnRoles else parsed

  role_assignments <- empty_role_assignments()
  known_roles <- names(role_assignments)

  for (role_name in known_roles) {
    role_value <- roles[[role_name]]
    if (!is.null(role_value)) {
      flattened <- as.character(unlist(role_value, use.names = FALSE))
      role_assignments[[role_name]] <- flattened[nzchar(flattened)]
    }
  }

  role_assignments
}

escape_regex_literal <- function(text_value) {
  gsub("([][{}()+*^$|\\?.])", "\\\\\\1", text_value)
}

build_measurement_groups_from_tissue_map <- function(tissue_map_candidates) {
  tissue_labels <- unname(unlist(tissue_map_candidates, use.names = FALSE))

  if (length(tissue_labels) == 0) {
    return(list(
      list(
        groupName = "defaultMeasurements",
        columnPrefixes = character(0),
        tier = 2
      )
    ))
  }

  lapply(tissue_labels, function(tissue_label) {
    list(
      groupName = paste0(tolower(tissue_label), "Cells"),
      columnPrefixes = c(paste0(tissue_label, "_")),
      tier = 2
    )
  })
}

build_generated_yaml_text <- function(generated_configuration) {
  guidance_lines <- c(
    "# Generated by build_minimal_configuration_from_schema().",
    "#",
    "# Overlap guidance:",
    "# - data.analysisTimepointColumn defaults to data.rawTimepointColumn when omitted.",
    "# - variables.exogenousVariables should include design variables such as group/timepoint.",
    "# - variables.stratifierVariables typically includes data.protectionColumn for outcome-driven strata.",
    "# - data.ignoredColumns captures known non-model fields to keep role assignment explicit."
  )

  paste(
    c(
      guidance_lines,
      "",
      yaml::as.yaml(generated_configuration, indent.mapping.sequence = TRUE)
    ),
    collapse = "\n"
  )
}

build_minimal_configuration_from_schema <- function(file_path, schema_summary = NULL, study_name = "networkR analysis", output_directory = "outputs") {
  if (is.null(schema_summary)) {
    schema_summary <- summarize_tabular_schema(file_path)
  }

  resolved_roles <- schema_summary$resolvedRoles
  generated_configuration <- list(
    study = list(
      studyName = study_name
    ),
    data = list(
      dataPath = dirname(file_path),
      filePattern = escape_regex_literal(basename(file_path)),
      inputFormat = schema_summary$suggestedInputFormat,
      idColumn = resolved_roles$idColumn,
      rawTimepointColumn = resolved_roles$rawTimepointColumn,
      analysisTimepointColumn = resolved_roles$analysisTimepointColumn,
      groupColumn = resolved_roles$groupColumn,
      protectionColumn = resolved_roles$protectionColumn,
      tissueColumn = resolved_roles$tissueColumn,
      ignoredColumns = resolved_roles$ignoredColumns,
      inferMeasurementColumns = TRUE,
      tissueMap = resolved_roles$tissueMapCandidates
    ),
    variables = list(
      exogenousVariables = resolved_roles$exogenousVariables,
      stratifierVariables = resolved_roles$stratifierVariables,
      derivedVariables = list(),
      measurementGroups = build_measurement_groups_from_tissue_map(resolved_roles$tissueMapCandidates)
    ),
    output = list(
      outputDirectory = output_directory
    )
  )

  list(
    configuration = generated_configuration,
    yamlText = build_generated_yaml_text(generated_configuration),
    unresolvedIssues = schema_summary$unresolvedIssues,
    reviewScaffold = schema_summary$reviewScaffold
  )
}