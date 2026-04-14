# networkR MCP server entry point.
#
# Typical client registrations:
#
# VS Code / GitHub Copilot / Positron:
# {
#   "servers": {
#     "networkR": {
#       "type": "stdio",
#       "command": "Rscript",
#       "args": ["/absolute/path/to/networkR/inst/mcp_server.R"]
#     }
#   }
# }
#
# Claude Desktop:
# {
#   "mcpServers": {
#     "networkR": {
#       "command": "Rscript",
#       "args": ["/absolute/path/to/networkR/inst/mcp_server.R"]
#     }
#   }
# }
#
# Claude Code:
# claude mcp add -s user networkR -- Rscript /absolute/path/to/networkR/inst/mcp_server.R
#
# Optional live-session registration:
# source("inst/mcp_server.R")
# networkR_mcp_session()

`%||%` <- function(left, right) {
  if (is.null(left)) right else left
}

absolute_path <- function(path) {
  grepl("^(/|[A-Za-z]:[/\\\\])", path)
}

script_file <- local({
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(NULL)
  }
  normalizePath(sub("^--file=", "", file_arg[[1]]), winslash = "/", mustWork = FALSE)
})

load_networkr_context <- function() {
  package_root <- NULL
  if (!is.null(script_file)) {
    candidate_root <- normalizePath(file.path(dirname(script_file), ".."), winslash = "/", mustWork = FALSE)
    if (file.exists(file.path(candidate_root, "DESCRIPTION")) && dir.exists(file.path(candidate_root, "R"))) {
      package_root <- candidate_root
    }
  }

  if (!is.null(package_root)) {
    package_env <- new.env(parent = baseenv())
    r_files <- sort(list.files(file.path(package_root, "R"), pattern = "[.]R$", full.names = TRUE))
    for (r_file in r_files) {
      sys.source(r_file, envir = package_env)
    }

    return(list(
      mode = "development",
      basePath = package_root,
      packageEnv = package_env,
      templatePath = normalizePath(file.path(package_root, "inst", "configuration_template.yml"), winslash = "/", mustWork = FALSE),
      examplePath = normalizePath(file.path(package_root, "inst", "example_configuration.yml"), winslash = "/", mustWork = FALSE)
    ))
  }

  if (!requireNamespace("networkR", quietly = TRUE)) {
    stop(
      paste(
        "Could not locate the networkR development checkout or an installed networkR package.",
        "Run this script from the repository, or install networkR before starting the MCP server."
      ),
      call. = FALSE
    )
  }

  list(
    mode = "installed",
    basePath = normalizePath(getwd(), winslash = "/", mustWork = FALSE),
    packageEnv = NULL,
    templatePath = system.file("configuration_template.yml", package = "networkR"),
    examplePath = system.file("example_configuration.yml", package = "networkR")
  )
}

.networkr_context <- load_networkr_context()
.networkr_state <- new.env(parent = emptyenv())
.networkr_state$last_pipeline_result <- NULL
.networkr_state$last_configuration_path <- NULL

require_mcp_packages <- function() {
  required_packages <- c("btw", "ellmer", "mcptools")
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      "Install required MCP packages before starting this server: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }
}

networkr_function <- function(name, exported = TRUE) {
  if (identical(.networkr_context$mode, "development")) {
    if (!exists(name, envir = .networkr_context$packageEnv, inherits = FALSE)) {
      stop("Could not locate function in development environment: ", name, call. = FALSE)
    }
    return(get(name, envir = .networkr_context$packageEnv, inherits = FALSE))
  }

  if (exported) {
    getExportedValue("networkR", name)
  } else {
    getFromNamespace(name, ns = "networkR")
  }
}

resolve_user_path <- function(path, base_path = .networkr_context$basePath) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("Path arguments must be a single non-empty string.", call. = FALSE)
  }

  expanded_path <- path.expand(path)
  if (!absolute_path(expanded_path)) {
    expanded_path <- file.path(base_path, expanded_path)
  }

  normalizePath(expanded_path, winslash = "/", mustWork = FALSE)
}

resolve_relative_to <- function(path, anchor_dir) {
  if (is.null(path) || !nzchar(path)) {
    return(path)
  }

  expanded_path <- path.expand(path)
  if (!absolute_path(expanded_path)) {
    expanded_path <- file.path(anchor_dir, expanded_path)
  }

  normalizePath(expanded_path, winslash = "/", mustWork = FALSE)
}

read_text_file <- function(path) {
  paste(suppressWarnings(readLines(path, warn = FALSE, encoding = "UTF-8")), collapse = "\n")
}

summarize_tabular_file <- function(file_path, max_unique_values) {
  tryCatch({
    summarize_tabular_schema <- networkr_function("summarize_tabular_schema", exported = FALSE)
    summarize_tabular_schema(file_path, max_unique_values = max_unique_values)
  }, error = function(error) {
    list(
      filePath = file_path,
      readSucceeded = FALSE,
      error = conditionMessage(error)
    )
  })
}

configuration_summary <- function(configuration, config_path = NULL) {
  measurement_groups <- configuration$variables$measurementGroups %||% list()
  stratifications <- configuration$stratifications %||% list()
  derived_variables <- configuration$variables$derivedVariables %||% list()

  list(
    configPath = config_path,
    studyName = configuration$study$studyName,
    backend = configuration$structureLearning$backend %||% "stagedtrees",
    dataPath = configuration$data$dataPath,
    dataPathExists = dir.exists(configuration$data$dataPath),
    outputDirectory = configuration$output$outputDirectory %||% NULL,
    exogenousVariables = configuration$variables$exogenousVariables,
    stratifierVariables = configuration$variables$stratifierVariables %||% character(0),
    derivedVariables = names(derived_variables),
    measurementGroups = lapply(measurement_groups, function(group_definition) {
      list(
        groupName = group_definition$groupName,
        tier = group_definition$tier,
        columnPrefixes = group_definition$columnPrefixes %||% character(0)
      )
    }),
    stratifications = lapply(stratifications, function(stratification_definition) {
      list(
        analysisName = stratification_definition$analysisName,
        variableName = stratification_definition$variableName,
        levels = stratification_definition$levels %||% character(0)
      )
    })
  )
}

load_configuration <- function(config_path) {
  resolved_config_path <- resolve_user_path(config_path)
  if (!file.exists(resolved_config_path)) {
    stop("Configuration file does not exist: ", resolved_config_path, call. = FALSE)
  }

  parse_configuration <- networkr_function("ParseConfiguration")
  validate_configuration <- networkr_function("ValidateConfiguration")

  configuration <- suppressWarnings(parse_configuration(resolved_config_path))
  config_dir <- dirname(resolved_config_path)
  configuration$data$dataPath <- resolve_relative_to(configuration$data$dataPath, config_dir)
  configuration$output$outputDirectory <- resolve_relative_to(configuration$output$outputDirectory %||% "outputs", config_dir)
  configuration <- validate_configuration(configuration)

  list(
    path = resolved_config_path,
    directory = config_dir,
    configuration = configuration
  )
}

count_averaged_edges <- function(structure_result, backend) {
  if (identical(backend, "stagedtrees")) {
    return(sum(structure_result$averagedNetwork != 0))
  }

  if (is.null(structure_result$averagedNetwork)) {
    return(0L)
  }

  nrow(as.data.frame(bnlearn::arcs(structure_result$averagedNetwork), stringsAsFactors = FALSE))
}

summarize_pipeline_result <- function(pipeline_result, config_path = NULL) {
  configuration <- pipeline_result$configuration
  structure_result <- pipeline_result$fullModelResult
  backend <- configuration$structureLearning$backend %||% "stagedtrees"
  output_directory <- configuration$output$outputDirectory %||% NULL
  saved_plot_files <- character(0)

  if (!is.null(output_directory) && dir.exists(output_directory)) {
    saved_plot_files <- list.files(
      output_directory,
      pattern = "[.](png|pdf|svg)$",
      full.names = TRUE,
      recursive = FALSE,
      ignore.case = TRUE
    )
  }

  list(
    configPath = config_path,
    studyName = configuration$study$studyName,
    backend = backend,
    rows = nrow(pipeline_result$subjectTable),
    measurementColumnCount = length(pipeline_result$imputationResult$measurementColumns),
    measurementColumns = utils::head(pipeline_result$imputationResult$measurementColumns, 12),
    activeExogenousVariables = structure_result$activeExogenousVariables,
    modelVariableCount = ncol(structure_result$modelData),
    bootstrapEdgeCount = nrow(structure_result$bootstrapStrength),
    averagedEdgeCount = count_averaged_edges(structure_result, backend),
    threshold = structure_result$threshold,
    primaryScore = if (length(structure_result$networkScores) > 0) as.numeric(structure_result$networkScores[[1]]) else NA_real_,
    stratifiedAnalyses = lapply(names(pipeline_result$plotObjects$stratifiedModels), function(analysis_name) {
      analysis_result <- pipeline_result$plotObjects$stratifiedModels[[analysis_name]]
      list(
        analysisName = analysis_name,
        levels = names(analysis_result$levelResults),
        comparisons = names(analysis_result$comparisons)
      )
    }),
    savePlots = isTRUE(configuration$output$savePlots),
    outputDirectory = output_directory,
    savedPlotFiles = saved_plot_files,
    storedResultName = "lastPipelineResult"
  )
}

networkr_get_configuration_schema <- function() {
  required_paths <- networkr_function("required_configuration_paths", exported = FALSE)()

  list(
    success = TRUE,
    mode = .networkr_context$mode,
    projectBasePath = .networkr_context$basePath,
    templatePath = .networkr_context$templatePath,
    templateYaml = read_text_file(.networkr_context$templatePath),
    exampleConfigurationPath = .networkr_context$examplePath,
    exampleConfigurationYaml = read_text_file(.networkr_context$examplePath),
    requiredPaths = vapply(required_paths, paste, collapse = ".", FUN.VALUE = character(1)),
    workflow = c(
      "Call networkr_explore_data_directory to inspect the input files — it reports suggestedInputFormat, grouped roleAssignments, unresolvedIssues, and a reviewScaffold per file.",
      "OR call networkr_generate_configuration_from_data to go directly from a single data file to a minimal YAML draft with a reviewScaffold embedded.",
      "Review the reviewScaffold: move column names between subjectId, sampleId, tissue, exogenousNotStratified, stratifiedBy, measurement, ignored, and unsure as needed.",
      "When the scaffold needs adjustments, call networkr_refine_configuration_from_scaffold with the edited columnRoles YAML block to regenerate the YAML from the corrected assignments.",
      "For 'wide_by_tissue': set data.tissueColumn and data.tissueMap. Omit data.measurementColumns when data.inferMeasurementColumns is true.",
      "Use networkr_write_configuration_file to save the YAML to disk.",
      "Use networkr_validate_configuration_file before running the pipeline.",
      "Run networkr_run_pipeline_from_config and then networkr_describe_last_pipeline_result for follow-up questions."
    )
  )
}

networkr_explore_data_directory <- function(data_path, max_unique_values = 10L) {
  tryCatch({
    resolved_data_path <- resolve_user_path(data_path)
    if (!dir.exists(resolved_data_path)) {
      return(list(success = FALSE, error = paste("Data directory does not exist:", resolved_data_path)))
    }

    data_files <- sort(list.files(
      resolved_data_path,
      pattern = "[.](csv|xls|xlsx)$",
      full.names = TRUE,
      recursive = FALSE,
      ignore.case = TRUE
    ))

    if (length(data_files) == 0L) {
      return(list(
        success = FALSE,
        dataPath = resolved_data_path,
        error = "No .csv, .xls, or .xlsx files were found in the requested directory."
      ))
    }

    list(
      success = TRUE,
      dataPath = resolved_data_path,
      files = lapply(data_files, summarize_tabular_file, max_unique_values = max_unique_values)
    )
  }, error = function(error) {
    list(success = FALSE, error = conditionMessage(error))
  })
}

networkr_validate_configuration_file <- function(config_path) {
  tryCatch({
    loaded_configuration <- load_configuration(config_path)
    list(
      success = TRUE,
      valid = TRUE,
      summary = configuration_summary(loaded_configuration$configuration, loaded_configuration$path)
    )
  }, error = function(error) {
    list(
      success = FALSE,
      valid = FALSE,
      configPath = resolve_user_path(config_path),
      error = conditionMessage(error)
    )
  })
}

networkr_generate_configuration_from_data <- function(data_file_path, study_name = "networkR analysis", output_directory = "outputs") {
  tryCatch({
    resolved_data_file <- resolve_user_path(data_file_path)
    if (!file.exists(resolved_data_file)) {
      return(list(
        success = FALSE,
        error = paste("Data file does not exist:", resolved_data_file)
      ))
    }

    summarize_tabular_schema <- networkr_function("summarize_tabular_schema", exported = FALSE)
    build_minimal_configuration_from_schema <- networkr_function("build_minimal_configuration_from_schema", exported = FALSE)
    build_configuration <- networkr_function("BuildConfiguration")

    schema_summary <- summarize_tabular_schema(resolved_data_file)
    generated <- build_minimal_configuration_from_schema(
      resolved_data_file,
      schema_summary = schema_summary,
      study_name = study_name,
      output_directory = output_directory
    )
    validated_configuration <- build_configuration(generated$configuration)

    list(
      success = TRUE,
      dataFilePath = resolved_data_file,
      summary = configuration_summary(validated_configuration),
      yamlText = generated$yamlText,
      reviewScaffold = generated$reviewScaffold,
      unresolvedIssues = generated$unresolvedIssues,
      roleAssignments = schema_summary$roleAssignments,
      resolvedRoles = schema_summary$resolvedRoles
    )
  }, error = function(error) {
    list(success = FALSE, error = conditionMessage(error))
  })
}

networkr_refine_configuration_from_scaffold <- function(scaffold_yaml, data_file_path, study_name = "networkR analysis", output_directory = "outputs") {
  tryCatch({
    resolved_data_file <- resolve_user_path(data_file_path)
    if (!file.exists(resolved_data_file)) {
      return(list(
        success = FALSE,
        error = paste("Data file does not exist:", resolved_data_file)
      ))
    }

    parse_role_scaffold_fn   <- networkr_function("parse_role_scaffold",   exported = FALSE)
    read_tabular_data_fn     <- networkr_function("read_tabular_data",     exported = FALSE)
    suggest_input_format_fn  <- networkr_function("suggest_input_format",  exported = FALSE)
    resolve_role_assignments_fn <- networkr_function("resolve_role_assignments", exported = FALSE)
    build_column_role_scaffold_fn <- networkr_function("build_column_role_scaffold", exported = FALSE)
    build_minimal_configuration_from_schema_fn <- networkr_function("build_minimal_configuration_from_schema", exported = FALSE)
    build_configuration <- networkr_function("BuildConfiguration")

    data_frame       <- read_tabular_data_fn(resolved_data_file)
    role_assignments <- parse_role_scaffold_fn(scaffold_yaml)
    resolved_roles   <- resolve_role_assignments_fn(data_frame, role_assignments)
    suggested_format <- suggest_input_format_fn(data_frame)

    schema_summary <- list(
      filePath           = resolved_data_file,
      nRows              = nrow(data_frame),
      nCols              = ncol(data_frame),
      suggestedInputFormat = suggested_format,
      roleAssignments    = role_assignments,
      resolvedRoles      = resolved_roles,
      unresolvedIssues   = list(),
      reviewScaffold     = build_column_role_scaffold_fn(list(
        columns         = lapply(colnames(data_frame), function(col_name) {
          list(
            name           = col_name,
            roleSuggestion = list(rationale = "User-provided assignment.")
          )
        }),
        roleAssignments = role_assignments
      ))
    )

    generated <- build_minimal_configuration_from_schema_fn(
      resolved_data_file,
      schema_summary   = schema_summary,
      study_name       = study_name,
      output_directory = output_directory
    )
    validated_configuration <- build_configuration(generated$configuration)

    list(
      success          = TRUE,
      dataFilePath     = resolved_data_file,
      yamlText         = generated$yamlText,
      reviewScaffold   = generated$reviewScaffold,
      unresolvedIssues = generated$unresolvedIssues,
      roleAssignments  = role_assignments,
      resolvedRoles    = resolved_roles,
      summary          = configuration_summary(validated_configuration)
    )
  }, error = function(error) {
    list(success = FALSE, error = conditionMessage(error))
  })
}

networkr_write_configuration_file <- function(config_path, yaml_text, overwrite = FALSE) {
  tryCatch({
    resolved_config_path <- resolve_user_path(config_path)
    if (file.exists(resolved_config_path) && !isTRUE(overwrite)) {
      return(list(
        success = FALSE,
        configPath = resolved_config_path,
        error = "Configuration file already exists. Re-run with overwrite = true to replace it."
      ))
    }

    dir.create(dirname(resolved_config_path), recursive = TRUE, showWarnings = FALSE)
    writeLines(enc2utf8(yaml_text), con = resolved_config_path, useBytes = TRUE)

    validation_result <- networkr_validate_configuration_file(resolved_config_path)
    list(
      success = TRUE,
      configPath = resolved_config_path,
      bytesWritten = unname(file.info(resolved_config_path)$size),
      validation = validation_result
    )
  }, error = function(error) {
    list(success = FALSE, error = conditionMessage(error))
  })
}

networkr_run_pipeline_from_config <- function(config_path, output_dir = NULL) {
  tryCatch({
    loaded_configuration <- load_configuration(config_path)
    configuration <- loaded_configuration$configuration

    if (!is.null(output_dir) && nzchar(output_dir)) {
      configuration$output$outputDirectory <- resolve_relative_to(output_dir, loaded_configuration$directory)
      configuration$output$savePlots <- TRUE
    }

    run_pipeline <- networkr_function("RunPipeline")
    pipeline_result <- run_pipeline(configuration)

    .networkr_state$last_pipeline_result <- pipeline_result
    .networkr_state$last_configuration_path <- loaded_configuration$path

    list(
      success = TRUE,
      runSucceeded = TRUE,
      summary = summarize_pipeline_result(pipeline_result, loaded_configuration$path)
    )
  }, error = function(error) {
    list(success = FALSE, runSucceeded = FALSE, error = conditionMessage(error))
  })
}

networkr_reload_package <- function() {
  if (!identical(.networkr_context$mode, "development")) {
    return(list(
      success = FALSE,
      error = "networkr_reload_package is only available in development mode (running from the package source directory)."
    ))
  }

  package_root <- .networkr_context$basePath
  r_files <- sort(list.files(file.path(package_root, "R"), pattern = "[.]R$", full.names = TRUE))

  loaded <- character(0)
  failed <- character(0)

  for (r_file in r_files) {
    tryCatch({
      sys.source(r_file, envir = .networkr_context$packageEnv)
      loaded <- c(loaded, basename(r_file))
    }, error = function(error) {
      failed <<- c(failed, paste0(basename(r_file), ": ", conditionMessage(error)))
    })
  }

  list(
    success = length(failed) == 0L,
    filesLoaded = loaded,
    filesFailed = failed,
    message = if (length(failed) == 0L) {
      paste("Reloaded", length(loaded), "R source files. All subsequent MCP calls will use the updated code.")
    } else {
      paste("Reloaded with errors.", length(failed), "files failed.")
    }
  )
}

networkr_describe_last_pipeline_result <- function() {
  if (is.null(.networkr_state$last_pipeline_result)) {
    return(list(
      success = FALSE,
      error = "No pipeline result is stored in this MCP server session yet. Run networkr_run_pipeline_from_config first."
    ))
  }

  list(
    success = TRUE,
    summary = summarize_pipeline_result(
      .networkr_state$last_pipeline_result,
      .networkr_state$last_configuration_path
    )
  )
}

networkR_mcp_tools <- function() {
  require_mcp_packages()

  c(
    btw::btw_tools("docs", "env", "sessioninfo", "pkg"),
    list(
      ellmer::tool(
        networkr_get_configuration_schema,
        name = "networkr_get_configuration_schema",
        description = paste(
          "Return the networkR configuration template, example configuration, required field paths,",
          "and workflow guidance for building YAML from a natural-language analysis goal.",
          "Use this before drafting a new configuration."
        ),
        arguments = list()
      ),
      ellmer::tool(
        networkr_explore_data_directory,
        name = "networkr_explore_data_directory",
        description = paste(
          "Inspect a directory of .csv, .xls, or .xlsx files and summarize each file's columns, classes,",
          "missingness, unique values, numeric ranges, suggested roles, unresolved issues, and a copy/paste review scaffold.",
          "Use this before mapping dataset columns into a networkR configuration."
        ),
        arguments = list(
          data_path = ellmer::type_string("Path to the directory that contains the input spreadsheets or CSV files for the analysis."),
          max_unique_values = ellmer::type_integer(
            description = "Maximum number of distinct non-missing example values to return for each column. Defaults to 10.",
            required = FALSE
          )
        )
      ),
      ellmer::tool(
        networkr_write_configuration_file,
        name = "networkr_write_configuration_file",
        description = paste(
          "Write YAML text to a configuration file on disk, optionally overwriting an existing file,",
          "then immediately validate the saved configuration. Use this when the model has produced a draft YAML file."
        ),
        arguments = list(
          config_path = ellmer::type_string("Destination path for the YAML configuration file, relative to the project root or absolute."),
          yaml_text = ellmer::type_string("The full YAML configuration text to write to disk."),
          overwrite = ellmer::type_boolean(
            description = "Set to true to replace an existing configuration file at the same path. Defaults to false.",
            required = FALSE
          )
        )
      ),
      ellmer::tool(
        networkr_generate_configuration_from_data,
        name = "networkr_generate_configuration_from_data",
        description = paste(
          "Inspect a single CSV or spreadsheet, infer column roles, and generate a minimal networkR YAML configuration draft.",
          "The response includes yamlText, a reviewScaffold, unresolved issues, and the resolved role assignments.",
          "Use this after reviewing networkr_explore_data_directory output when you want a compact starting config."
        ),
        arguments = list(
          data_file_path = ellmer::type_string("Path to the input CSV, XLS, or XLSX file to convert into a minimal networkR configuration draft."),
          study_name = ellmer::type_string(
            description = "Study name to embed in the generated YAML. Defaults to 'networkR analysis'.",
            required = FALSE
          ),
          output_directory = ellmer::type_string(
            description = "Output directory to place into the generated YAML. Defaults to 'outputs'.",
            required = FALSE
          )
        )
      ),
      ellmer::tool(
        networkr_refine_configuration_from_scaffold,
        name = "networkr_refine_configuration_from_scaffold",
        description = paste(
          "Accept an edited columnRoles YAML scaffold block, re-resolve the column assignments against the original data file,",
          "and regenerate the networkR YAML configuration. Use this when the user has corrected a reviewScaffold",
          "returned by networkr_generate_configuration_from_data or networkr_explore_data_directory.",
          "The scaffold_yaml argument should contain the full columnRoles: block with corrected column placements."
        ),
        arguments = list(
          scaffold_yaml = ellmer::type_string(
            "The edited columnRoles YAML block produced by reviewing a reviewScaffold. Must contain a 'columnRoles:' key at the top level."
          ),
          data_file_path = ellmer::type_string("Path to the original CSV, XLS, or XLSX data file used to generate the scaffold."),
          study_name = ellmer::type_string(
            description = "Study name to embed in the regenerated YAML. Defaults to 'networkR analysis'.",
            required = FALSE
          ),
          output_directory = ellmer::type_string(
            description = "Output directory to embed in the regenerated YAML. Defaults to 'outputs'.",
            required = FALSE
          )
        )
      ),
      ellmer::tool(
        networkr_validate_configuration_file,
        name = "networkr_validate_configuration_file",
        description = paste(
          "Parse a networkR YAML configuration file, normalize relative paths, validate required fields,",
          "and return a concise summary or a validation error. Use this before running the pipeline."
        ),
        arguments = list(
          config_path = ellmer::type_string("Path to an existing networkR YAML configuration file.")
        )
      ),
      ellmer::tool(
        networkr_run_pipeline_from_config,
        name = "networkr_run_pipeline_from_config",
        description = paste(
          "Run the full networkR pipeline from a validated YAML configuration file.",
          "Optionally override the output directory; when provided, plots are saved there.",
          "The most recent result is stored in the MCP server session for later inspection."
        ),
        arguments = list(
          config_path = ellmer::type_string("Path to the networkR YAML configuration file to execute."),
          output_dir = ellmer::type_string(
            description = "Optional output directory for saved plots and artifacts. Relative paths are resolved from the configuration file directory.",
            required = FALSE
          )
        )
      ),
      ellmer::tool(
        networkr_describe_last_pipeline_result,
        name = "networkr_describe_last_pipeline_result",
        description = paste(
          "Describe the most recent pipeline result stored in the MCP server session, including backend, sample count,",
          "measurement columns, averaged edge count, stratified analyses, and saved plot files.",
          "Use this for follow-up questions after running the pipeline."
        ),
        arguments = list()
      ),
      ellmer::tool(
        networkr_reload_package,
        name = "networkr_reload_package",
        description = paste(
          "Reload all R source files from the networkR development checkout into the MCP server's cached package environment.",
          "Use this after editing any R/ files to pick up code changes without restarting VS Code.",
          "Only available in development mode (running from the package source directory)."
        ),
        arguments = list()
      )
    )
  )
}

networkR_mcp_session <- function() {
  require_mcp_packages()
  btw::btw_mcp_session()
}

networkR_mcp_server <- function() {
  require_mcp_packages()
  btw::btw_mcp_server(tools = networkR_mcp_tools())
}

if (sys.nframe() == 0L) {
  networkR_mcp_server()
}