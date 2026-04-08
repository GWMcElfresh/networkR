collapse_bootstrap_strength <- function(bootstrap_strength, label_name) {
  bootstrap_strength |>
    dplyr::mutate(edge = paste(pmin(.data$from, .data$to), pmax(.data$from, .data$to), sep = " -- ")) |>
    dplyr::group_by(.data$edge) |>
    dplyr::summarise(strength = max(.data$strength), .groups = "drop") |>
    dplyr::mutate(label = label_name)
}

signed_bootstrap_strength <- function(bootstrap_strength, continuous_data, label_name) {
  collapsed_strength <- collapse_bootstrap_strength(bootstrap_strength, label_name)

  correlation_values <- vapply(collapsed_strength$edge, function(edge_name) {
    edge_nodes <- strsplit(edge_name, " -- ", fixed = TRUE)[[1]]
    if (!all(edge_nodes %in% colnames(continuous_data))) {
      return(NA_real_)
    }

    first_values <- suppressWarnings(as.numeric(continuous_data[[edge_nodes[[1]]]]))
    second_values <- suppressWarnings(as.numeric(continuous_data[[edge_nodes[[2]]]]))
    complete_rows <- stats::complete.cases(data.frame(first_values, second_values))
    if (sum(complete_rows) < 5) {
      return(NA_real_)
    }

    stats::cor(first_values[complete_rows], second_values[complete_rows], method = "pearson")
  }, numeric(1))

  collapsed_strength |>
    dplyr::mutate(
      correlation = correlation_values,
      signedStrength = ifelse(is.na(.data$correlation), .data$strength, sign(.data$correlation) * .data$strength)
    )
}

#' Compare two learned Bayesian network structures
#'
#' @param firstStructure A result from `LearnBayesianNetwork()`.
#' @param secondStructure A result from `LearnBayesianNetwork()`.
#' @param firstLabel Label for the first structure.
#' @param secondLabel Label for the second structure.
#'
#' @return A list containing a comparison table.
CompareStructures <- function(firstStructure, secondStructure, firstLabel = "First", secondLabel = "Second") {
  first_strength <- signed_bootstrap_strength(
    firstStructure$bootstrapStrength,
    firstStructure$continuousData,
    firstLabel
  )
  second_strength <- signed_bootstrap_strength(
    secondStructure$bootstrapStrength,
    secondStructure$continuousData,
    secondLabel
  )

  measurement_groups <- c(firstStructure$measurementGroups, secondStructure$measurementGroups)
  active_exogenous_variables <- unique(c(
    firstStructure$activeExogenousVariables,
    secondStructure$activeExogenousVariables
  ))

  comparison_table <- dplyr::full_join(
    first_strength |>
      dplyr::select(.data$edge, firstStrength = .data$strength, firstSignedStrength = .data$signedStrength),
    second_strength |>
      dplyr::select(.data$edge, secondStrength = .data$strength, secondSignedStrength = .data$signedStrength),
    by = "edge"
  ) |>
    tidyr::replace_na(list(firstStrength = 0, secondStrength = 0, firstSignedStrength = 0, secondSignedStrength = 0)) |>
    dplyr::mutate(
      delta = .data$firstSignedStrength - .data$secondSignedStrength,
      absoluteDelta = abs(.data$delta),
      direction = dplyr::case_when(
        .data$delta > 0.1 ~ paste0("Stronger in ", firstLabel),
        .data$delta < -0.1 ~ paste0("Stronger in ", secondLabel),
        TRUE ~ "Similar"
      ),
      edgeType = vapply(.data$edge, function(edge_name) {
        edge_nodes <- strsplit(edge_name, " -- ", fixed = TRUE)[[1]]
        classify_edge_type(edge_nodes[[1]], edge_nodes[[2]], active_exogenous_variables, measurement_groups)
      }, character(1))
    ) |>
    dplyr::arrange(dplyr::desc(.data$absoluteDelta))

  list(
    comparisonTable = comparison_table,
    firstLabel = firstLabel,
    secondLabel = secondLabel
  )
}