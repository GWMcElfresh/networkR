run_stratified_models <- function(discretization_result, configuration) {
  stratified_results <- list()

  for (stratification_definition in configuration$stratifications) {
    variable_name <- stratification_definition$variableName
    level_names <- stratification_definition$levels
    filter_definitions <- stratification_definition$filters %||% list()
    analysis_name <- stratification_definition$analysisName
    eligible_rows <- apply_filters(discretization_result$imputedData, filter_definitions)

    level_results <- list()
    level_plots <- list()
    level_event_tree_plots <- list()
    level_stage_split_heatmaps <- list()
    for (level_name in level_names) {
      level_rows <- eligible_rows[discretization_result$imputedData[[variable_name]][eligible_rows] == level_name]
      if (length(level_rows) == 0) {
        next
      }

      level_discretization_result <- subset_discretization_result(discretization_result, level_rows)
      level_structure_result <- LearnBayesianNetwork(level_discretization_result, configuration)
      level_results[[level_name]] <- level_structure_result
      level_plots[[level_name]] <- PlotBayesianNetwork(
        level_structure_result,
        configuration,
        title = paste(configuration$study$studyName, analysis_name, level_name, sep = " - ")
      )
      if (!is.null(level_structure_result$sevtModel)) {
        level_event_tree_plots[[level_name]] <- tryCatch(
          PlotEventTree(
            level_structure_result,
            configuration,
            title = paste(configuration$study$studyName, analysis_name, level_name, "event tree", sep = " - ")
          ),
          error = function(e) {
            message("PlotEventTree() failed for '", level_name, "': ", conditionMessage(e))
            NULL
          }
        )
        level_stage_split_heatmaps[[level_name]] <- tryCatch(
          PlotStageSplitHeatmap(
            level_structure_result,
            configuration,
            title = paste(configuration$study$studyName, analysis_name, level_name, "stage splits", sep = " - ")
          ),
          error = function(e) {
            message("PlotStageSplitHeatmap() failed for '", level_name, "': ", conditionMessage(e))
            NULL
          }
        )
      }
    }

    comparison_results <- list()
    comparison_plots <- list()
    if (length(level_results) >= 2) {
      level_pairs <- utils::combn(names(level_results), 2, simplify = FALSE)
      for (level_pair in level_pairs) {
        comparison_name <- paste(level_pair[[1]], level_pair[[2]], sep = "_vs_")
        comparison_result <- CompareStructures(
          level_results[[level_pair[[1]]]],
          level_results[[level_pair[[2]]]],
          firstLabel = level_pair[[1]],
          secondLabel = level_pair[[2]]
        )
        comparison_results[[comparison_name]] <- comparison_result
        comparison_plots[[comparison_name]] <- build_comparison_plot(
          comparison_result,
          configuration,
          paste(configuration$study$studyName, analysis_name, comparison_name, sep = " - ")
        )
      }
    }

    stratified_results[[analysis_name]] <- list(
      levelResults         = level_results,
      levelPlots           = level_plots,
      levelEventTreePlots  = level_event_tree_plots,
      levelStageSplitHeatmaps = level_stage_split_heatmaps,
      comparisons          = comparison_results,
      comparisonPlots      = comparison_plots
    )
  }

  stratified_results
}

save_pipeline_plots <- function(plot_objects, configuration) {
  output_directory <- configuration$output$outputDirectory
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)

  if (inherits(plot_objects$fullModel, "ggplot")) {
    SavePlot(plot_objects$fullModel, file.path(output_directory, "full_model_network.png"), configuration)
  }

  if (inherits(plot_objects$fullModelEventTree, "ggplot")) {
    SavePlot(plot_objects$fullModelEventTree, file.path(output_directory, "full_model_event_tree.png"), configuration)
  }

  if (inherits(plot_objects$fullModelStageSplitHeatmap, "ggplot")) {
    SavePlot(plot_objects$fullModelStageSplitHeatmap, file.path(output_directory, "full_model_stage_split_heatmap.png"), configuration)
  }

  for (analysis_name in names(plot_objects$stratifiedModels)) {
    analysis_plots <- plot_objects$stratifiedModels[[analysis_name]]
    for (level_name in names(analysis_plots$levelPlots)) {
      level_plot <- analysis_plots$levelPlots[[level_name]]
      if (inherits(level_plot, "ggplot")) {
        SavePlot(
          level_plot,
          file.path(output_directory, paste0(make_file_name_fragment(analysis_name), "_", make_file_name_fragment(level_name), "_network.png")),
          configuration
        )
      }
    }
    for (level_name in names(analysis_plots$levelEventTreePlots)) {
      tree_plot <- analysis_plots$levelEventTreePlots[[level_name]]
      if (inherits(tree_plot, "ggplot")) {
        SavePlot(
          tree_plot,
          file.path(output_directory, paste0(make_file_name_fragment(analysis_name), "_", make_file_name_fragment(level_name), "_event_tree.png")),
          configuration
        )
      }
    }
    for (level_name in names(analysis_plots$levelStageSplitHeatmaps)) {
      heatmap_plot <- analysis_plots$levelStageSplitHeatmaps[[level_name]]
      if (inherits(heatmap_plot, "ggplot")) {
        SavePlot(
          heatmap_plot,
          file.path(output_directory, paste0(make_file_name_fragment(analysis_name), "_", make_file_name_fragment(level_name), "_stage_split_heatmap.png")),
          configuration
        )
      }
    }
    for (comparison_name in names(analysis_plots$comparisonPlots)) {
      comparison_plot <- analysis_plots$comparisonPlots[[comparison_name]]
      if (inherits(comparison_plot, "ggplot")) {
        SavePlot(
          comparison_plot,
          file.path(output_directory, paste0(make_file_name_fragment(analysis_name), "_", make_file_name_fragment(comparison_name), "_comparison.png")),
          configuration
        )
      }
    }
  }
}

#' Run the analysis pipeline
#'
#' @param configuration A configuration object or YAML file path.
#'
#' @return A pipeline result object.
RunPipeline <- function(configuration) {
  if (is.character(configuration) && length(configuration) == 1) {
    configuration <- ParseConfiguration(configuration)
  } else {
    ValidateConfiguration(configuration)
  }

  subject_table <- ReadPanels(configuration)
  subject_table <- AddDerivedVariables(subject_table, configuration)
  imputation_result <- ImputeMeasurements(subject_table, configuration)
  discretization_result <- DiscretizeMeasurements(imputation_result, configuration)
  full_model_result <- LearnBayesianNetwork(discretization_result, configuration)

  full_model_network_plot <- PlotBayesianNetwork(full_model_result, configuration, paste0(configuration$study$studyName, " - Full model"))
  full_model_event_tree_plot <- if (!is.null(full_model_result$sevtModel)) {
    tryCatch(
      PlotEventTree(full_model_result, configuration, paste0(configuration$study$studyName, " - Full model event tree")),
      error = function(e) {
        message("PlotEventTree() failed for full model: ", conditionMessage(e))
        NULL
      }
    )
  } else {
    NULL
  }
  full_model_stage_split_heatmap <- if (!is.null(full_model_result$sevtModel)) {
    tryCatch(
      PlotStageSplitHeatmap(full_model_result, configuration, paste0(configuration$study$studyName, " - Full model stage splits")),
      error = function(e) {
        message("PlotStageSplitHeatmap() failed for full model: ", conditionMessage(e))
        NULL
      }
    )
  } else {
    NULL
  }

  plot_objects <- list(
    fullModel          = full_model_network_plot,
    fullModelEventTree = full_model_event_tree_plot,
    fullModelStageSplitHeatmap = full_model_stage_split_heatmap,
    stratifiedModels   = run_stratified_models(discretization_result, configuration)
  )

  pipeline_result <- structure(
    list(
      configuration = configuration,
      subjectTable = subject_table,
      imputationResult = imputation_result,
      discretizationResult = discretization_result,
      fullModelResult = full_model_result,
      plotObjects = plot_objects
    ),
    class = "networkRPipelineResult"
  )

  if (isTRUE(configuration$output$savePlots)) {
    save_pipeline_plots(plot_objects, configuration)
  }

  pipeline_result
}

print.networkRPipelineResult <- function(x, ...) {
  backend <- x$configuration$structureLearning$backend %||% "stagedtrees"
  cores   <- x$configuration$structureLearning$bootstrapCores %||% 1
  cat("networkR pipeline result\n")
  cat("  Study name:        ", x$configuration$study$studyName, "\n", sep = "")
  cat("  Backend:           ", backend, "\n", sep = "")
  cat("  Bootstrap cores:   ", cores, "\n", sep = "")
  cat("  Rows:              ", nrow(x$subjectTable), "\n", sep = "")
  cat("  Measurement cols:  ", length(x$imputationResult$measurementColumns), "\n", sep = "")
  cat("  Stratified models: ", length(x$plotObjects$stratifiedModels), "\n", sep = "")
  invisible(x)
}

summary.networkRPipelineResult <- function(object, ...) {
  print.networkRPipelineResult(object, ...)
}