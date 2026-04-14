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

# For the stagedtrees backend, extract edges from the averaged adjacency matrix
# (which is a plain matrix) and attach bootstrap arc frequencies as strength.
build_plot_tables_stagedtrees <- function(structure_result) {
  averaged_adj <- structure_result$averagedNetwork  # matrix
  arc_indices <- which(averaged_adj > 0, arr.ind = TRUE)

  if (nrow(arc_indices) == 0) {
    edge_table <- data.frame(from = character(0), to = character(0), strength = numeric(0),
                             association = numeric(0), stringsAsFactors = FALSE)
    node_names <- rownames(averaged_adj) %||% colnames(averaged_adj)
  } else {
    bs_lookup <- build_edge_strength_lookup(structure_result$bootstrapStrength)
    assoc_lookup <- compute_edge_association(structure_result$sevtModel)

    edge_table <- data.frame(
      from = rownames(averaged_adj)[arc_indices[, 1]],
      to   = colnames(averaged_adj)[arc_indices[, 2]],
      stringsAsFactors = FALSE
    )
    edge_table$strength <- vapply(seq_len(nrow(edge_table)), function(i) {
      key <- paste(edge_table$from[[i]], edge_table$to[[i]], sep = " -> ")
      v <- bs_lookup[key]; if (is.na(v)) 0 else v
    }, numeric(1))
    edge_table$association <- vapply(seq_len(nrow(edge_table)), function(i) {
      key <- paste(edge_table$from[[i]], edge_table$to[[i]], sep = " -> ")
      v <- assoc_lookup[key]; if (length(v) == 0 || is.na(v)) 0 else v
    }, numeric(1))

    node_names <- unique(c(edge_table$from, edge_table$to))
  }
  list(edge_table = edge_table, node_names = node_names)
}

build_plot_tables_bnlearn <- function(structure_result) {
  edge_table <- as.data.frame(bnlearn::arcs(structure_result$averagedNetwork), stringsAsFactors = FALSE)
  bs_lookup  <- build_edge_strength_lookup(structure_result$bootstrapStrength)

  if (nrow(edge_table) == 0) {
    node_names <- bnlearn::nodes(structure_result$averagedNetwork)
    edge_table <- data.frame(from = character(0), to = character(0),
                             strength = numeric(0), association = numeric(0))
  } else {
    node_names <- unique(c(edge_table$from, edge_table$to))
    edge_table$strength <- vapply(seq_len(nrow(edge_table)), function(i) {
      key <- paste(edge_table$from[[i]], edge_table$to[[i]], sep = " -> ")
      v <- bs_lookup[key]; if (is.na(v)) 0 else v
    }, numeric(1))
    edge_table$association <- edge_table$strength  # use bootstrap strength as association proxy
  }
  list(edge_table = edge_table, node_names = node_names)
}

# Compute per-edge total variation distance between conditional distributions
# across all stage pairs of the parent variable in a fitted sevt model.
# Returns a named vector: names are "from -> to", values are mean TV in [0,1].
compute_edge_association <- function(sevt_model) {
  if (is.null(sevt_model) || !inherits(sevt_model, "sevt")) {
    return(character(0))
  }

  parent_list <- tryCatch(
    suppressWarnings(suppressMessages(stagedtrees::as_parentslist(sevt_model, silent = TRUE))),
    error = function(e) NULL
  )
  if (is.null(parent_list)) return(character(0))

  result <- list()

  for (child_var in names(parent_list)) {
    entry   <- parent_list[[child_var]]
    parents <- entry$parents
    if (is.null(parents) || length(parents) == 0) next

    child_probs <- sevt_model$prob[[child_var]]
    if (is.null(child_probs)) next
    stages <- sevt_model$stages[[child_var]]
    if (is.null(stages)) next

    unobserved <- sevt_model$name_unobserved %||% character(0)
    unique_stages <- unique(stages[!stages %in% c(unobserved, NA)])
    if (length(unique_stages) < 2) next

    # For each parent, mark edge association as mean TV over all stage pairs
    for (parent_var in parents) {
      stage_probs <- lapply(unique_stages, function(s) {
        p <- tryCatch(child_probs[[s]], error = function(e) NULL)
        if (is.null(p)) return(NULL)
        as.numeric(p)
      })
      stage_probs <- Filter(Negate(is.null), stage_probs)
      if (length(stage_probs) < 2) next

      pairs <- utils::combn(length(stage_probs), 2, simplify = FALSE)
      tv_values <- vapply(pairs, function(pair) {
        p1 <- stage_probs[[pair[[1]]]]
        p2 <- stage_probs[[pair[[2]]]]
        if (length(p1) != length(p2)) return(0)
        0.5 * sum(abs(p1 - p2))
      }, numeric(1))

      key <- paste(parent_var, child_var, sep = " -> ")
      result[[key]] <- mean(tv_values, na.rm = TRUE)
    }
  }

  unlist(result)
}

build_stage_split_summary <- function(structure_result) {
  if (is.null(structure_result$sevtModel)) {
    stop("build_stage_split_summary() requires a stagedtrees structure result.", call. = FALSE)
  }

  stage_info <- summary(structure_result$sevtModel)$stages.info
  variable_order <- structure_result$variableOrder %||% names(stage_info)

  summary_rows <- lapply(seq_along(variable_order), function(variable_index) {
    variable_name <- variable_order[[variable_index]]
    stage_table <- stage_info[[variable_name]]

    if (is.null(stage_table)) {
      return(NULL)
    }

    stage_table <- as.data.frame(stage_table, stringsAsFactors = FALSE)
    probability_columns <- setdiff(colnames(stage_table), c("stage", "npaths", "sample.size"))
    if (length(probability_columns) == 0) {
      return(NULL)
    }

    stage_table$variable <- variable_name
    stage_table$variableIndex <- variable_index
    stage_table$stage <- as.character(stage_table$stage)
    stage_table$stageRank <- seq_len(nrow(stage_table))
    stage_table$npaths <- as.numeric(stage_table$npaths)
    stage_table$sampleSize <- as.numeric(stage_table[["sample.size"]] %||% rep(NA_real_, nrow(stage_table)))
    stage_table$positionWeight <- ifelse(
      is.na(stage_table$sampleSize),
      NA_real_,
      stage_table$sampleSize / pmax(stage_table$npaths, 1)
    )
    stage_table$nodeGroup <- node_group_name(
      variable_name,
      structure_result$activeExogenousVariables,
      structure_result$measurementGroups
    )

    tidyr::pivot_longer(
      stage_table,
      cols = dplyr::all_of(probability_columns),
      names_to = "split",
      values_to = "probability"
    )
  })

  dplyr::bind_rows(summary_rows)
}

build_position_stage_map <- function(ceg_object) {
  variable_names <- names(ceg_object$tree)
  position_stage_map <- list()

  for (variable_name in variable_names) {
    unique_positions <- unique(ceg_object$positions[[variable_name]])
    for (position_name in unique_positions) {
      stage_index <- which.max(ceg_object$positions[[variable_name]] == position_name)
      stage_name <- ceg_object$stages[[variable_name]][stage_index]
      position_stage_map[[paste0(variable_name, ":", position_name)]] <- stage_name
    }
  }

  position_stage_map
}

build_stage_metric_lookup <- function(stage_summary, metric_name) {
  metric_table <- stage_summary |>
    dplyr::distinct(.data$variable, .data$stage, .data[[metric_name]]) |>
    dplyr::mutate(key = paste(.data$variable, .data$stage, sep = "::"))

  stats::setNames(metric_table[[metric_name]], metric_table$key)
}

build_plot_tables <- function(structure_result, configuration) {
  is_stagedtrees_backend <- !is.null(structure_result$sevtModel)

  if (is_stagedtrees_backend) {
    tables <- build_plot_tables_stagedtrees(structure_result)
  } else {
    tables <- build_plot_tables_bnlearn(structure_result)
  }

  edge_table <- tables$edge_table
  node_names <- tables$node_names

  node_table <- data.frame(name = node_names, stringsAsFactors = FALSE)
  node_table$nodeGroup <- vapply(
    node_table$name, node_group_name, character(1),
    structure_result$activeExogenousVariables, structure_result$measurementGroups
  )
  node_colours <- lapply(node_table$nodeGroup, resolve_group_colours,
    active_exogenous_variables = structure_result$activeExogenousVariables,
    measurement_groups = structure_result$measurementGroups,
    plotting_configuration = configuration$plotting
  )
  node_table$fillColour   <- vapply(node_colours, `[[`, character(1), "fillColour")
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

#' Plot an averaged network graph with group-colored nodes and strength-weighted edges
#'
#' Produces a ggraph DAG plot. Nodes are colored by measurement group (and exogenous status).
#' Edge width and alpha encode arc association strength (bootstrap frequency for bnlearn;
#' mean total-variation distance for stagedtrees).
#'
#' @param structureResult A result from `LearnBayesianNetwork()`.
#' @param configuration A networkR configuration object.
#' @param title Optional plot title.
#'
#' @return A ggplot object.
PlotBayesianNetwork <- function(structureResult, configuration, title = NULL) {
  plot_tables <- build_plot_tables(structureResult, configuration)

  # Use association column if present, otherwise fall back to strength
  edge_df <- plot_tables$edgeTable
  if (!"association" %in% colnames(edge_df)) {
    edge_df$association <- edge_df$strength
  }

  graph_object <- igraph::graph_from_data_frame(
    edge_df,
    directed = TRUE,
    vertices = plot_tables$nodeTable
  )

  title_text <- title %||% paste0(configuration$study$studyName, " - averaged network")

  ggraph::ggraph(graph_object, layout = configuration$plotting$networkLayout) +
    ggraph::geom_edge_link(
      ggplot2::aes(edge_width = .data$association, edge_alpha = .data$association),
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
    ggraph::scale_edge_width(range = c(0.3, 2.5), name = "Association strength") +
    ggraph::scale_edge_alpha(range = c(0.2, 0.85), guide = "none") +
    ggplot2::labs(title = title_text) +
    egg::theme_article(base_size = configuration$plotting$baseSize) +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}

#' Plot stage-to-split probabilities for a staged event tree
#'
#' Produces a heatmap where fill encodes conditional split probability and point
#' size encodes the stage sample size.
#'
#' @param structureResult A result from `LearnBayesianNetwork()` with a `$sevtModel` slot.
#' @param configuration A networkR configuration object.
#' @param title Optional plot title.
#'
#' @return A ggplot object.
PlotStageSplitHeatmap <- function(structureResult, configuration, title = NULL) {
  if (is.null(structureResult$sevtModel)) {
    stop("PlotStageSplitHeatmap() requires a stagedtrees backend result with a $sevtModel slot.", call. = FALSE)
  }

  stage_summary <- build_stage_split_summary(structureResult)
  if (nrow(stage_summary) == 0) {
    return(NULL)
  }

  stage_levels <- stage_summary |>
    dplyr::distinct(.data$variable, .data$variableIndex, .data$stage, .data$stageRank, .data$sampleSize) |>
    dplyr::arrange(.data$variableIndex, .data$stageRank) |>
    dplyr::mutate(
      stageDisplay = paste0(
        "Stage ", .data$stage,
        " (n=", ifelse(is.na(.data$sampleSize), "NA", format(round(.data$sampleSize, 1), trim = TRUE)), ")"
      )
    )

  stage_summary <- dplyr::left_join(
    stage_summary,
    stage_levels |>
      dplyr::select("variable", "stage", "stageDisplay"),
    by = c("variable", "stage")
  )

  stage_summary$variable <- factor(
    stage_summary$variable,
    levels = unique(stage_levels$variable)
  )
  stage_summary$stageDisplay <- factor(
    stage_summary$stageDisplay,
    levels = rev(unique(stage_levels$stageDisplay))
  )

  title_text <- title %||% paste0(configuration$study$studyName, " - stage to split probabilities")

  plot_object <- ggplot2::ggplot(
    stage_summary,
    ggplot2::aes(x = .data$split, y = .data$stageDisplay)
  ) +
    ggplot2::geom_tile(
      ggplot2::aes(fill = .data$probability),
      colour = "white",
      linewidth = 0.3
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = sprintf("%.2f", .data$probability)),
      size = 3
    ) +
    ggplot2::facet_grid(variable ~ ., scales = "free_y", space = "free_y") +
    ggplot2::scale_fill_gradient(
      low = "grey95",
      high = "steelblue4",
      limits = c(0, 1),
      name = "P(split | stage)"
    ) +
    ggplot2::labs(
      title = title_text,
      x = "Split",
      y = NULL
    ) +
    egg::theme_article(base_size = configuration$plotting$baseSize) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1)
    )

  if (any(!is.na(stage_summary$sampleSize))) {
    plot_object <- plot_object +
      ggplot2::geom_point(
        ggplot2::aes(size = .data$sampleSize),
        shape = 21,
        fill = NA,
        colour = "black",
        stroke = 0.3,
        alpha = 0.7
      ) +
      ggplot2::scale_size_area(max_size = 7, name = "Stage sample size")
  }

  plot_object
}

# Convert a value from make_stages_col (integer palette index or character
# colour string) to a character colour string.
palette_or_colour <- function(x) {
  if (is.character(x)) return(x)
  pal <- grDevices::palette()
  pal[((as.integer(x) - 1L) %% length(pal)) + 1L]
}

# Build an igraph from a stagedtrees ceg object without using ceg2adjmat
# (which is not exported in all stagedtrees versions). Edge labels record the
# split value(s) that produce each transition.
build_ceg_graph_components <- function(ceg_obj) {
  variable_names <- names(ceg_obj$tree)
  pos_labeled <- lapply(seq_along(ceg_obj$positions), function(index) {
    paste0(variable_names[[index]], ":", ceg_obj$positions[[variable_names[[index]]]])
  })
  names(pos_labeled) <- variable_names

  edge_rows <- list()
  edge_index <- 1L

  if (length(variable_names) >= 2) {
    for (variable_index in seq(2, length(variable_names))) {
      parent_var <- variable_names[[variable_index - 1L]]
      child_var <- variable_names[[variable_index]]
      branch_values <- as.character(ceg_obj$tree[[parent_var]])
      n_levels <- length(branch_values)
      parent_positions <- unique(pos_labeled[[parent_var]])

      for (parent_position in parent_positions) {
        first_match <- min(which(pos_labeled[[parent_var]] == parent_position))
        child_range <- ((first_match - 1L) * n_levels + 1L):(first_match * n_levels)
        child_positions <- pos_labeled[[child_var]][child_range]

        for (branch_index in seq_along(child_positions)) {
          edge_rows[[edge_index]] <- data.frame(
            from = parent_position,
            to = child_positions[[branch_index]],
            parentVariable = parent_var,
            childVariable = child_var,
            edgeValue = branch_values[[branch_index]],
            stringsAsFactors = FALSE
          )
          edge_index <- edge_index + 1L
        }
      }
    }
  }

  if (length(variable_names) > 0) {
    last_var <- variable_names[[length(variable_names)]]
    branch_values <- as.character(ceg_obj$tree[[last_var]])

    for (parent_position in unique(pos_labeled[[last_var]])) {
      for (branch_value in branch_values) {
        edge_rows[[edge_index]] <- data.frame(
          from = parent_position,
          to = "END",
          parentVariable = last_var,
          childVariable = "END",
          edgeValue = branch_value,
          stringsAsFactors = FALSE
        )
        edge_index <- edge_index + 1L
      }
    }
  }

  if (length(edge_rows) == 0) {
    edge_table <- data.frame(
      from = character(0),
      to = character(0),
      parentVariable = character(0),
      childVariable = character(0),
      edgeWeight = numeric(0),
      edgeLabel = character(0),
      stringsAsFactors = FALSE
    )
  } else {
    edge_table <- dplyr::bind_rows(edge_rows) |>
      dplyr::group_by(.data$from, .data$to, .data$parentVariable, .data$childVariable) |>
      dplyr::summarise(
        edgeWeight = dplyr::n(),
        edgeLabel = paste(unique(.data$edgeValue), collapse = " / "),
        .groups = "drop"
      )
  }

  vertex_table <- data.frame(
    name = c(unique(unlist(pos_labeled, use.names = FALSE)), "END"),
    stringsAsFactors = FALSE
  )

  graph_object <- igraph::graph_from_data_frame(
    edge_table,
    directed = TRUE,
    vertices = vertex_table
  )

  list(
    graph = graph_object,
    edgeTable = edge_table,
    positions = pos_labeled
  )
}

build_event_tree_axis_table <- function(layout_data) {
  as.data.frame(layout_data) |>
    dplyr::mutate(variable = sub(":.*$", "", .data$name)) |>
    dplyr::filter(.data$variable != "END") |>
    dplyr::group_by(.data$variable) |>
    dplyr::summarise(y = mean(.data$y), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(.data$y))
}

build_event_tree_edge_label_table <- function(layout_data, graph_object) {
  edge_table <- igraph::as_data_frame(graph_object, what = "edges")
  if (nrow(edge_table) == 0) {
    return(data.frame(
      x = numeric(0),
      y = numeric(0),
      edgeLabel = character(0),
      parentVariable = character(0),
      stringsAsFactors = FALSE
    ))
  }

  layout_table <- as.data.frame(layout_data)
  from_index <- match(edge_table$from, layout_table$name)
  to_index <- match(edge_table$to, layout_table$name)

  data.frame(
    x = (layout_table$x[from_index] + layout_table$x[to_index]) / 2,
    y = (layout_table$y[from_index] + layout_table$y[to_index]) / 2,
    edgeLabel = edge_table$edgeLabel,
    parentVariable = edge_table$parentVariable,
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(nzchar(.data$edgeLabel))
}

#' Plot a staged event tree (Chain Event Graph) with stage-colored nodes
#'
#' Produces a ggraph visualization of the CEG derived from a fitted stagedtrees
#' model. Nodes are colored by stage identity (stagedtrees convention), and node
#' borders are colored by measurement group (networkR convention). Requires
#' `configuration$structureLearning$backend == "stagedtrees"`.
#'
#' @param structureResult A result from `LearnBayesianNetwork()` with a `$sevtModel`
#'   slot (i.e. from the `stagedtrees` backend).
#' @param configuration A networkR configuration object.
#' @param title Optional plot title.
#' @param colorBy One of `"stage"` (color nodes by stage; default) or `"group"`
#'   (color nodes by measurement group - same as `PlotBayesianNetwork`).
#'
#' @return A ggplot object.
PlotEventTree <- function(structureResult, configuration, title = NULL, colorBy = "stage") {
  if (is.null(structureResult$sevtModel)) {
    stop("PlotEventTree() requires a stagedtrees backend result with a $sevtModel slot.", call. = FALSE)
  }

  sevt_model <- structureResult$sevtModel
  title_text  <- title %||% paste0(configuration$study$studyName, " - event tree")

  # Guard: ceg() enumerates every path through the tree, so models with many
  # variables (or large branching factors) produce astronomically large position
  # vectors.  Refuse gracefully above a safe threshold.
  n_tree_vars   <- length(sevt_model$tree)
  total_paths   <- prod(vapply(sevt_model$tree, length, integer(1)))
  path_limit    <- configuration$plotting$cegPathLimit %||% 50000L
  if (total_paths > path_limit) {
    message(
      "PlotEventTree(): skipping - the staged event tree has ", n_tree_vars,
      " variables with ", total_paths, " paths (limit: ", path_limit, "). ",
      "Consider setting a smaller model or raising configuration$plotting$cegPathLimit."
    )
    return(NULL)
  }

  # Build CEG object then convert to igraph via adjacency matrix.
  # We inline the ceg -> adjacency-matrix conversion to avoid depending on
  # stagedtrees:::ceg2adjmat, which is not exported in all package versions.
  ceg_obj     <- stagedtrees::ceg(sevt_model)
  graph_components <- build_ceg_graph_components(ceg_obj)
  ceg_igraph  <- graph_components$graph
  stage_summary <- build_stage_split_summary(structureResult)
  position_stage_map <- build_position_stage_map(ceg_obj)
  stage_weight_lookup <- build_stage_metric_lookup(stage_summary, "positionWeight")
  stage_size_lookup <- build_stage_metric_lookup(stage_summary, "sampleSize")

  # Node names are "<variable>:<position>" e.g. "Gender:1", plus "END"
  node_names <- igraph::V(ceg_igraph)$name

  # Split node names into variable and position parts
  node_variable <- sub(":.*$", "", node_names)

  # --- stage fill colours ---
  if (identical(colorBy, "stage")) {
    stage_col_list <- stagedtrees::make_stages_col(ceg_obj)

    stage_fill <- vapply(node_names, function(nm) {
      if (nm == "END") return("grey60")
      stg <- position_stage_map[[nm]]
      var <- sub(":.*$", "", nm)
      col_var <- stage_col_list[[var]]
      if (!is.null(col_var) && !is.null(stg) && stg %in% names(col_var)) {
        palette_or_colour(col_var[[stg]])
      } else {
        "grey80"
      }
    }, character(1))

  } else {
    # colorBy == "group": fill by measurement group of the variable
    stage_fill <- vapply(node_variable, function(vname) {
      if (vname == "END") return("grey60")
      grp <- node_group_name(vname,
                             structureResult$activeExogenousVariables,
                             structureResult$measurementGroups)
      resolve_group_colours(grp, structureResult$activeExogenousVariables,
                            structureResult$measurementGroups,
                            configuration$plotting)$fillColour
    }, character(1))
  }

  node_weight <- vapply(node_names, function(node_name) {
    if (node_name == "END") {
      return(0)
    }

    variable_name <- sub(":.*$", "", node_name)
    stage_name <- position_stage_map[[node_name]]
    lookup_key <- paste(variable_name, stage_name, sep = "::")
    metric_value <- stage_weight_lookup[lookup_key]
    if (length(metric_value) == 0 || is.na(metric_value)) 0 else unname(metric_value)
  }, numeric(1))
  node_sample_size <- vapply(node_names, function(node_name) {
    if (node_name == "END") {
      return(NA_real_)
    }

    variable_name <- sub(":.*$", "", node_name)
    stage_name <- position_stage_map[[node_name]]
    lookup_key <- paste(variable_name, stage_name, sep = "::")
    metric_value <- stage_size_lookup[lookup_key]
    if (length(metric_value) == 0 || is.na(metric_value)) NA_real_ else unname(metric_value)
  }, numeric(1))
  frequency_labels <- ifelse(
    is.na(node_sample_size) | node_sample_size <= 0,
    "",
    paste0("n=", format(round(node_sample_size, 1), trim = TRUE))
  )
  unobserved_stage_names <- unique(c(sevt_model$name_unobserved %||% character(0), "UNOBSERVED"))
  stage_labels <- vapply(node_names, function(node_name) {
    if (node_name == "END") {
      return("END")
    }

    stage_name <- position_stage_map[[node_name]]
    if (length(stage_name) == 0 || is.na(stage_name)) {
      return("Start")
    }

    if (stage_name %in% unobserved_stage_names) {
      return("U")
    }

    paste0("S", stage_name)
  }, character(1))

  # Border always encodes measurement group
  group_border <- vapply(node_variable, function(vname) {
    if (vname == "END") return("grey40")
    grp <- node_group_name(vname,
                           structureResult$activeExogenousVariables,
                           structureResult$measurementGroups)
    resolve_group_colours(grp, structureResult$activeExogenousVariables,
                          structureResult$measurementGroups,
                          configuration$plotting)$borderColour
  }, character(1))

  igraph::V(ceg_igraph)$stageFill   <- unname(stage_fill)
  igraph::V(ceg_igraph)$groupBorder <- unname(group_border)
  igraph::V(ceg_igraph)$label       <- stage_labels
  igraph::V(ceg_igraph)$nodeWeight  <- node_weight
  igraph::V(ceg_igraph)$sampleSize  <- node_sample_size
  igraph::V(ceg_igraph)$frequencyLabel <- frequency_labels

  edge_ends <- igraph::ends(ceg_igraph, igraph::E(ceg_igraph))
  edge_density <- node_weight[match(edge_ends[, 2], node_names)]
  edge_density[is.na(edge_density)] <- 0
  igraph::E(ceg_igraph)$transitionDensity <- edge_density

  layout_data <- ggraph::create_layout(ceg_igraph, layout = "tree")
  axis_table <- build_event_tree_axis_table(layout_data)
  edge_label_table <- build_event_tree_edge_label_table(layout_data, ceg_igraph)

  show_frequency <- isTRUE(configuration$plotting$eventTreeShowFrequency %||% TRUE) && any(node_weight > 0, na.rm = TRUE)

  if (show_frequency) {
    plot_object <- ggraph::ggraph(layout_data) +
      ggraph::geom_edge_link(
        ggplot2::aes(edge_width = .data$transitionDensity, edge_alpha = .data$transitionDensity),
        arrow = grid::arrow(type = "open", length = grid::unit(2.5, "mm")),
        end_cap = ggraph::circle(4, "mm"),
        colour = "grey40",
        show.legend = TRUE
      ) +
      ggraph::geom_node_point(
        ggplot2::aes(
          fill = .data$stageFill,
          colour = .data$groupBorder,
          size = .data$nodeWeight,
          alpha = .data$nodeWeight
        ),
        shape = 21,
        stroke = 1.2,
        show.legend = TRUE
      ) +
      ggraph::geom_node_text(
        ggplot2::aes(label = .data$label),
        vjust = -1.0,
        size = 2.4
      ) +
      ggraph::geom_node_text(
        ggplot2::aes(label = .data$frequencyLabel),
        vjust = 1.5,
        size = 2.1,
        colour = "grey25"
      )
  } else {
    plot_object <- ggraph::ggraph(layout_data) +
      ggraph::geom_edge_link(
        arrow = grid::arrow(type = "open", length = grid::unit(2.5, "mm")),
        end_cap = ggraph::circle(4, "mm"),
        colour = "grey40",
        width  = 0.6
      ) +
      ggraph::geom_node_point(
        ggplot2::aes(fill = .data$stageFill, colour = .data$groupBorder),
        shape = 21,
        size = 6,
        stroke = 1.2,
        show.legend = FALSE
      ) +
      ggraph::geom_node_text(
        ggplot2::aes(label = .data$label),
        size = 2.4
      )
  }

  if (nrow(edge_label_table) > 0) {
    plot_object <- plot_object +
      ggrepel::geom_label_repel(
        data = edge_label_table,
        ggplot2::aes(x = .data$x, y = .data$y, label = .data$edgeLabel),
        inherit.aes = FALSE,
        fill = grDevices::adjustcolor("white", alpha.f = 0.8),
        colour = "grey15",
        size = 2.6,
        label.padding = grid::unit(0.12, "lines"),
        label.r = grid::unit(0.08, "lines"),
        box.padding = grid::unit(0.1, "lines"),
        point.padding = grid::unit(0.05, "lines"),
        min.segment.length = 0,
        segment.size = 0.2,
        segment.alpha = 0.4,
        max.overlaps = Inf,
        seed = configuration$study$randomSeed
      )
  }

  plot_object +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_colour_identity() +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.05, 0.18))) +
    ggplot2::scale_y_continuous(
      breaks = axis_table$y,
      labels = axis_table$variable,
      expand = ggplot2::expansion(mult = c(0.03, 0.08))
    ) +
    ggplot2::labs(title = title_text, x = NULL, y = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    egg::theme_article(base_size = configuration$plotting$baseSize) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(face = "bold", colour = "grey20"),
      plot.margin = grid::unit(c(8, 28, 8, 12), "pt")
    ) +
    {
      if (show_frequency) ggplot2::scale_size_continuous(range = c(3, 10), name = "Position density") else NULL
    } +
    {
      if (show_frequency) ggplot2::scale_alpha_continuous(range = c(0.45, 1), guide = "none") else NULL
    } +
    {
      if (show_frequency) ggraph::scale_edge_width(range = c(0.3, 2.4), name = "Transition density") else NULL
    } +
    {
      if (show_frequency) ggraph::scale_edge_alpha(range = c(0.2, 0.8), guide = "none") else NULL
    }
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
