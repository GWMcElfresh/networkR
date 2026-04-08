resolve_group_colours <- function(node_group, active_exogenous_variables, measurement_groups, plotting_configuration) {
  if (identical(node_group, "exogenous")) {
    return(plotting_configuration$groupColours$exogenous)
  }

  for (measurement_group in measurement_groups) {
    if (identical(node_group, measurement_group$groupName)) {
      return(list(fillColour = measurement_group$fillColour, borderColour = measurement_group$borderColour))
    }
  }

  plotting_configuration$groupColours$defaultMeasurement
}

build_edge_strength_lookup <- function(bootstrap_strength) {
  edge_table <- bootstrap_strength |>
    dplyr::mutate(edgeKey = paste(.data$from, .data$to, sep = " -> "))
  stats::setNames(edge_table$strength, edge_table$edgeKey)
}

build_plot_tables <- function(structure_result, configuration) {
  edge_table <- as.data.frame(bnlearn::arcs(structure_result$averagedNetwork), stringsAsFactors = FALSE)
  edge_strength_lookup <- build_edge_strength_lookup(structure_result$bootstrapStrength)

  if (nrow(edge_table) == 0) {
    node_names <- bnlearn::nodes(structure_result$averagedNetwork)
    edge_table <- data.frame(from = character(0), to = character(0), strength = numeric(0))
  } else {
    node_names <- unique(c(edge_table$from, edge_table$to))
    edge_table$strength <- vapply(seq_len(nrow(edge_table)), function(row_index) {
      edge_key <- paste(edge_table$from[[row_index]], edge_table$to[[row_index]], sep = " -> ")
      edge_strength_lookup[[edge_key]] %||% 0
    }, numeric(1))
  }

  node_table <- data.frame(name = node_names, stringsAsFactors = FALSE)
  node_table$nodeGroup <- vapply(node_table$name, node_group_name, character(1), structure_result$activeExogenousVariables, structure_result$measurementGroups)
  node_colours <- lapply(node_table$nodeGroup, resolve_group_colours,
    active_exogenous_variables = structure_result$activeExogenousVariables,
    measurement_groups = structure_result$measurementGroups,
    plotting_configuration = configuration$plotting
  )
  node_table$fillColour <- vapply(node_colours, `[[`, character(1), "fillColour")
  node_table$borderColour <- vapply(node_colours, `[[`, character(1), "borderColour")

  list(edgeTable = edge_table, nodeTable = node_table)
}

build_comparison_plot <- function(comparison_result, configuration, title_text) {
  comparison_table <- comparison_result$comparisonTable

  plot_object <- ggplot2::ggplot(
    comparison_table,
    ggplot2::aes(
      x = .data$secondSignedStrength,
      y = .data$firstSignedStrength,
      colour = .data$direction,
      label = .data$edge
    )
  ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", colour = "grey40") +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey40") +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_point(ggplot2::aes(size = .data$absoluteDelta), alpha = 0.75) +
    ggrepel::geom_label_repel(
      data = comparison_table[comparison_table$absoluteDelta > configuration$plotting$edgeLabelThreshold, , drop = FALSE],
      size = 3.5,
      max.overlaps = 10,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = c(
      stats::setNames("steelblue", paste0("Stronger in ", comparison_result$firstLabel)),
      stats::setNames("firebrick", paste0("Stronger in ", comparison_result$secondLabel)),
      Similar = "grey50"
    )) +
    ggplot2::scale_size_continuous(range = c(1, 5)) +
    ggplot2::labs(
      title = title_text,
      x = paste0("Signed arc strength - ", comparison_result$secondLabel),
      y = paste0("Signed arc strength - ", comparison_result$firstLabel),
      colour = NULL,
      size = "|delta|"
    ) +
    egg::theme_article(base_size = configuration$plotting$baseSize)

  if ("edgeType" %in% colnames(comparison_table)) {
    plot_object <- plot_object + ggplot2::facet_grid(edgeType ~ .)
  }

  plot_object
}

#' Plot an averaged Bayesian network with ggraph
#'
#' @param structureResult A result from `LearnBayesianNetwork()`.
#' @param configuration A networkR configuration object.
#' @param title Optional plot title.
#'
#' @return A ggplot object.
PlotBayesianNetwork <- function(structureResult, configuration, title = NULL) {
  plot_tables <- build_plot_tables(structureResult, configuration)
  graph_object <- igraph::graph_from_data_frame(
    plot_tables$edgeTable,
    directed = TRUE,
    vertices = plot_tables$nodeTable
  )

  title_text <- title %||% paste0(
    configuration$study$studyName,
    " - averaged network"
  )

  ggraph::ggraph(graph_object, layout = configuration$plotting$networkLayout) +
    ggraph::geom_edge_link(
      ggplot2::aes(edge_width = .data$strength, edge_alpha = .data$strength),
      arrow = grid::arrow(type = "closed", length = grid::unit(3, "mm")),
      end_cap = ggraph::circle(3, "mm"),
      show.legend = TRUE
    ) +
    ggraph::geom_node_point(
      ggplot2::aes(fill = .data$fillColour, colour = .data$borderColour),
      shape = 21,
      size = 7,
      stroke = 1.1,
      show.legend = FALSE
    ) +
    ggraph::geom_node_text(
      ggplot2::aes(label = .data$name),
      repel = TRUE,
      size = 3.5
    ) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_colour_identity() +
    ggraph::scale_edge_width(range = c(0.3, 2.5), name = "Bootstrap strength") +
    ggraph::scale_edge_alpha(range = c(0.2, 0.85), guide = "none") +
    ggplot2::labs(title = title_text) +
    egg::theme_article(base_size = configuration$plotting$baseSize) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' Save a plot using package defaults
#'
#' @param plotObject A ggplot object.
#' @param filePath Output file path.
#' @param configuration Optional networkR configuration object.
#'
#' @return The saved file path, invisibly.
SavePlot <- function(plotObject, filePath, configuration = NULL) {
  if (is.null(configuration)) {
    configuration <- BuildConfiguration()
  }
  save_configuration <- configuration$plotting$save
  ggplot2::ggsave(
    filename = filePath,
    plot = plotObject,
    width = save_configuration$width,
    height = save_configuration$height,
    dpi = save_configuration$dpi,
    bg = save_configuration$bg,
    device = save_configuration$device
  )

  invisible(filePath)
}