default_configuration <- function() {
  list(
    study = list(
      studyName = "networkR analysis",
      randomSeed = 42
    ),
    data = list(
      dataPath = "data_sheets",
      filePattern = "[.](xlsx|csv)$",
      idColumn = "SubjectId",
      rawTimepointColumn = "TimepointLabel2",
      analysisTimepointColumn = "Timepoint",
      groupColumn = "GroupName",
      protectionColumn = "Protection",
      tissueColumn = "Tissue",
      cellTypeColumn = "CellType",
      fractionColumn = "Fraction",
      cellTypePrefix = "Pct_",
      tissueMap = list(
        liver = "Liver",
        pbmc = "PBMC",
        bm = "BM"
      ),
      cellTypeRenames = list(
        list(
          match = "Promyelocytes",
          replace = "Neutrophils",
          exceptTissue = "BM"
        )
      )
    ),
    variables = list(
      exogenousVariables = c("GroupName", "VaccinationStatus", "Timepoint"),
      stratifierVariables = c("Protection"),
      derivedVariables = list(
        VaccinationStatus = list(
          type = "ifEquals",
          sourceColumn = "GroupName",
          equalsValue = "Control Animal",
          trueValue = "Unvaccinated",
          falseValue = "Vaccinated"
        ),
        Protection = list(
          type = "replaceTargetWhenEquals",
          sourceColumn = "GroupName",
          targetColumn = "Protection",
          equalsValue = "Control Animal",
          replacementValue = NA_character_
        )
      ),
      measurementGroups = list(
        list(
          groupName = "liverCells",
          columnPrefixes = c("Liver_"),
          tier = 2,
          fillColour = "lightgreen",
          borderColour = "darkgreen"
        ),
        list(
          groupName = "pbmcCells",
          columnPrefixes = c("PBMC_"),
          tier = 3,
          fillColour = "lightblue",
          borderColour = "steelblue"
        ),
        list(
          groupName = "boneMarrowCells",
          columnPrefixes = c("BM_"),
          tier = 3,
          fillColour = "navajowhite",
          borderColour = "sienna4"
        )
      )
    ),
    imputation = list(
      method = "pmm",
      numberOfImputations = 5,
      maximumIterations = 20,
      completeDataAction = 1,
      excludedPredictors = c("SubjectId", "Protection")
    ),
    discretization = list(
      method = "hartemink",
      initialBreaks = 20,
      finalBreaks = 3,
      fallbackMethod = "quantile",
      fallbackBreaks = 2
    ),
    structureLearning = list(
      algorithms = c("hillClimbing", "tabu"),
      score = "bic",
      selectedBootstrapAlgorithm = "tabu",
      bootstrapReplicates = 1000,
      averageNetworkThreshold = 0.4,
      parameterMethod = "bayes",
      imaginarySampleSize = 1
    ),
    stratifications = list(),
    plotting = list(
      themeName = "article",
      baseSize = 13,
      networkLayout = "fr",
      edgeLabelThreshold = 0.2,
      comparisonBarThreshold = 0.15,
      groupColours = list(
        exogenous = list(
          fillColour = "mistyrose",
          borderColour = "firebrick"
        ),
        defaultMeasurement = list(
          fillColour = "grey85",
          borderColour = "grey30"
        )
      ),
      save = list(
        width = 13,
        height = 8,
        dpi = 300,
        bg = "white",
        device = "png"
      )
    ),
    output = list(
      savePlots = FALSE,
      outputDirectory = "outputs"
    )
  )
}

required_configuration_paths <- function() {
  list(
    c("study", "studyName"),
    c("study", "randomSeed"),
    c("data", "dataPath"),
    c("data", "idColumn"),
    c("data", "rawTimepointColumn"),
    c("data", "analysisTimepointColumn"),
    c("data", "groupColumn"),
    c("data", "protectionColumn"),
    c("data", "cellTypeColumn"),
    c("data", "fractionColumn"),
    c("variables", "exogenousVariables"),
    c("variables", "measurementGroups"),
    c("imputation", "method"),
    c("discretization", "method"),
    c("structureLearning", "algorithms"),
    c("plotting", "networkLayout")
  )
}

extract_nested_value <- function(configuration, path_components) {
  current_value <- configuration
  for (path_component in path_components) {
    if (!is.list(current_value) || is.null(current_value[[path_component]])) {
      return(NULL)
    }
    current_value <- current_value[[path_component]]
  }
  current_value
}

normalize_configuration <- function(configuration) {
  class(configuration) <- unique(c("networkRConfiguration", class(configuration)))
  configuration
}

#' Build a networkR configuration object
#'
#' @param overrides A named list of configuration overrides.
#'
#' @return A validated configuration object.
BuildConfiguration <- function(overrides = list()) {
  if (!is.list(overrides)) {
    stop("overrides must be a list.", call. = FALSE)
  }

  configuration <- utils::modifyList(default_configuration(), overrides)

  # utils::modifyList merges named list elements recursively but silently ignores
  # unnamed list elements (YAML sequences / R list(list(...)) patterns).
  # Re-apply fields that use unnamed lists directly when provided in overrides.
  if ("stratifications" %in% names(overrides)) {
    configuration$stratifications <- overrides$stratifications
  }
  if (!is.null(overrides[["variables"]])) {
    if ("measurementGroups" %in% names(overrides$variables)) {
      configuration$variables$measurementGroups <- overrides$variables$measurementGroups
    }
    if ("derivedVariables" %in% names(overrides$variables)) {
      configuration$variables$derivedVariables <- overrides$variables$derivedVariables
    }
  }
  if (!is.null(overrides[["data"]])) {
    if ("cellTypeRenames" %in% names(overrides$data)) {
      configuration$data$cellTypeRenames <- overrides$data$cellTypeRenames
    }
  }

  configuration <- normalize_configuration(configuration)
  ValidateConfiguration(configuration)
  configuration
}

#' Parse a YAML configuration file
#'
#' @param filePath Path to a YAML configuration file.
#'
#' @return A validated configuration object.
ParseConfiguration <- function(filePath) {
  if (!file.exists(filePath)) {
    stop("Configuration file does not exist: ", filePath, call. = FALSE)
  }

  parsed_configuration <- yaml::read_yaml(filePath)
  if (is.null(parsed_configuration)) {
    parsed_configuration <- list()
  }

  BuildConfiguration(parsed_configuration)
}

#' Validate a networkR configuration object
#'
#' @param configuration A configuration list.
#'
#' @return The validated configuration object, invisibly.
ValidateConfiguration <- function(configuration) {
  for (required_path in required_configuration_paths()) {
    required_value <- extract_nested_value(configuration, required_path)
    if (is.null(required_value)) {
      stop(
        "Missing required configuration field: ",
        paste(required_path, collapse = "."),
        call. = FALSE
      )
    }
  }

  if (!is.numeric(configuration$study$randomSeed) || length(configuration$study$randomSeed) != 1) {
    stop("study.randomSeed must be a single numeric value.", call. = FALSE)
  }

  if (!is.character(configuration$variables$exogenousVariables) || length(configuration$variables$exogenousVariables) == 0) {
    stop("variables.exogenousVariables must be a non-empty character vector.", call. = FALSE)
  }

  if (!is.list(configuration$variables$measurementGroups)) {
    stop("variables.measurementGroups must be a list.", call. = FALSE)
  }

  for (measurement_group in configuration$variables$measurementGroups) {
    if (is.null(measurement_group$groupName) || is.null(measurement_group$tier)) {
      stop("Each measurement group must include groupName and tier.", call. = FALSE)
    }
  }

  invisible(normalize_configuration(configuration))
}

print.networkRConfiguration <- function(x, ...) {
  cat("networkR configuration\n")
  cat("  Study name: ", x$study$studyName, "\n", sep = "")
  cat("  Data path: ", x$data$dataPath, "\n", sep = "")
  cat(
    "  Exogenous variables: ",
    paste(x$variables$exogenousVariables, collapse = ", "),
    "\n",
    sep = ""
  )
  cat(
    "  Measurement groups: ",
    paste(vapply(x$variables$measurementGroups, `[[`, character(1), "groupName"), collapse = ", "),
    "\n",
    sep = ""
  )
  invisible(x)
}