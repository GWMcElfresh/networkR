read_panel_file <- function(file_path, configuration) {
  file_extension <- tolower(tools::file_ext(file_path))

  if (identical(file_extension, "csv")) {
    panel_table <- utils::read.csv(file_path, check.names = FALSE, stringsAsFactors = FALSE)
  } else {
    panel_table <- as.data.frame(readxl::read_excel(file_path), check.names = FALSE)
  }

  valid_name_mask <- nchar(trimws(colnames(panel_table))) > 0 &
    !grepl("^\\.\\.\\.\\d+$", colnames(panel_table))
  panel_table <- panel_table[, valid_name_mask, drop = FALSE]

  required_columns <- c(
    configuration$data$idColumn,
    configuration$data$rawTimepointColumn,
    configuration$data$groupColumn,
    configuration$data$protectionColumn,
    configuration$data$cellTypeColumn,
    configuration$data$fractionColumn
  )
  assert_has_columns(panel_table, required_columns, object_name = basename(file_path))

  tissue_name <- parse_tissue_name(file_path, configuration$data$tissueMap)

  if (configuration$data$rawTimepointColumn != configuration$data$analysisTimepointColumn) {
    colnames(panel_table)[colnames(panel_table) == configuration$data$rawTimepointColumn] <- configuration$data$analysisTimepointColumn
  }

  panel_table[[configuration$data$idColumn]] <- as.character(panel_table[[configuration$data$idColumn]])
  cleaned_cell_types <- clean_cell_type(
    panel_table[[configuration$data$cellTypeColumn]],
    configuration$data$cellTypePrefix
  )
  prefixed_cell_types <- paste0(tissue_name, "_", cleaned_cell_types)
  prefixed_cell_types <- apply_cell_type_renames(
    prefixed_cell_types,
    tissue_name,
    configuration$data$cellTypeRenames,
    configuration$data$cellTypePrefix
  )

  panel_table[[configuration$data$cellTypeColumn]] <- prefixed_cell_types
  panel_table[[configuration$data$tissueColumn]] <- tissue_name

  identifier_columns <- c(
    configuration$data$idColumn,
    configuration$data$analysisTimepointColumn,
    configuration$data$tissueColumn,
    configuration$data$protectionColumn,
    configuration$data$groupColumn,
    configuration$data$cellTypeColumn
  )
  fraction_column <- configuration$data$fractionColumn

  summarized_table <- panel_table |>
    dplyr::group_by(dplyr::across(dplyr::all_of(identifier_columns))) |>
    dplyr::summarise(
      measurementValue = sum(.data[[fraction_column]], na.rm = TRUE),
      .groups = "drop"
    )

  summarized_table |>
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(c(
        configuration$data$idColumn,
        configuration$data$analysisTimepointColumn,
        configuration$data$tissueColumn,
        configuration$data$protectionColumn,
        configuration$data$groupColumn
      )),
      names_from = dplyr::all_of(configuration$data$cellTypeColumn),
      values_from = "measurementValue"
    )
}

#' Read panel files into a joined subject table
#'
#' @param configuration A networkR configuration object.
#'
#' @return A joined subject-level table.
ReadPanels <- function(configuration) {
  ValidateConfiguration(configuration)

  panel_files <- list.files(
    configuration$data$dataPath,
    pattern = configuration$data$filePattern,
    full.names = TRUE,
    recursive = FALSE
  )

  if (length(panel_files) == 0) {
    stop("No panel files found in dataPath: ", configuration$data$dataPath, call. = FALSE)
  }

  panel_tables <- lapply(panel_files, read_panel_file, configuration = configuration)
  join_keys <- c(
    configuration$data$idColumn,
    configuration$data$analysisTimepointColumn,
    configuration$data$tissueColumn,
    configuration$data$protectionColumn,
    configuration$data$groupColumn
  )

  Reduce(function(left_table, right_table) {
    dplyr::full_join(left_table, right_table, by = join_keys)
  }, panel_tables)
}

apply_derived_variable <- function(subject_table, variable_name, rule_definition) {
  if (identical(rule_definition$type, "ifEquals")) {
    source_values <- subject_table[[rule_definition$sourceColumn]]
    subject_table[[variable_name]] <- ifelse(
      source_values == rule_definition$equalsValue,
      rule_definition$trueValue,
      rule_definition$falseValue
    )
    return(subject_table)
  }

  if (identical(rule_definition$type, "replaceTargetWhenEquals")) {
    source_values <- subject_table[[rule_definition$sourceColumn]]
    target_values <- as.character(subject_table[[rule_definition$targetColumn]])
    replacement_value <- rule_definition$replacementValue
    replacement_value <- if (is.null(replacement_value)) NA_character_ else replacement_value
    target_values[source_values == rule_definition$equalsValue] <- replacement_value
    subject_table[[variable_name]] <- target_values
    return(subject_table)
  }

  stop("Unsupported derived variable rule type: ", rule_definition$type, call. = FALSE)
}

#' Add derived design variables to a subject table
#'
#' @param subjectTable A subject-level table.
#' @param configuration A networkR configuration object.
#'
#' @return The subject table with derived variables added.
AddDerivedVariables <- function(subjectTable, configuration) {
  subject_table <- subjectTable
  derived_variables <- configuration$variables$derivedVariables %||% list()

  for (variable_name in names(derived_variables)) {
    subject_table <- apply_derived_variable(subject_table, variable_name, derived_variables[[variable_name]])
  }

  subject_table
}