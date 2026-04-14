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

    first_values  <- suppressWarnings(as.numeric(continuous_data[[edge_nodes[[1]]]]))
    second_values <- suppressWarnings(as.numeric(continuous_data[[edge_nodes[[2]]]]))
    complete_rows <- stats::complete.cases(data.frame(first_values, second_values))
    if (sum(complete_rows) < 5) {
      return(NA_real_)
    }

    stats::cor(first_values[complete_rows], second_values[complete_rows], method = "pearson")
  }, numeric(1))

  collapsed_strength |>
    dplyr::mutate(
      correlation    = correlation_values,
      signedStrength = ifelse(is.na(.data$correlation), .data$strength, sign(.data$correlation) * .data$strength)
    )
}

normalize_directed_bootstrap_strength <- function(bootstrap_strength, label_name) {
  required_columns <- c("from", "to", "strength")
  missing_columns <- setdiff(required_columns, colnames(bootstrap_strength))
  if (length(missing_columns) > 0) {
    stop(
      "bootstrapStrength is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  bootstrap_strength |>
    dplyr::transmute(
      edge = paste(.data$from, .data$to, sep = " -> "),
      strength = as.numeric(.data$strength),
      label = label_name
    ) |>
    dplyr::group_by(.data$edge, .data$label) |>
    dplyr::summarise(strength = max(.data$strength), .groups = "drop")
}

compare_directed_bootstrap_strength <- function(first_bootstrap_strength, second_bootstrap_strength, first_label, second_label, tolerance) {
  first_table <- normalize_directed_bootstrap_strength(first_bootstrap_strength, first_label)
  second_table <- normalize_directed_bootstrap_strength(second_bootstrap_strength, second_label)

  edge_table <- dplyr::full_join(
    first_table |>
      dplyr::select("edge", firstStrength = "strength"),
    second_table |>
      dplyr::select("edge", secondStrength = "strength"),
    by = "edge"
  ) |>
    tidyr::replace_na(list(firstStrength = 0, secondStrength = 0)) |>
    dplyr::mutate(
      absoluteDifference = abs(.data$firstStrength - .data$secondStrength),
      exceededTolerance = .data$absoluteDifference > tolerance,
      presentInFirst = .data$firstStrength > 0,
      presentInSecond = .data$secondStrength > 0,
      presentInBoth = .data$presentInFirst & .data$presentInSecond
    ) |>
    dplyr::arrange(dplyr::desc(.data$absoluteDifference), .data$edge)

  edge_count <- nrow(edge_table)
  exceedance_count <- sum(edge_table$exceededTolerance)
  shared_nonzero_count <- sum(edge_table$presentInBoth)

  list(
    edgeTable = edge_table,
    summary = list(
      edgeCount = edge_count,
      meanAbsoluteDifference = if (edge_count == 0) 0 else mean(edge_table$absoluteDifference),
      maxAbsoluteDifference = if (edge_count == 0) 0 else max(edge_table$absoluteDifference),
      exceedanceCount = exceedance_count,
      exceedanceRate = if (edge_count == 0) 0 else exceedance_count / edge_count,
      sharedNonzeroEdgeCount = shared_nonzero_count,
      sharedNonzeroEdgeRate = if (edge_count == 0) 0 else shared_nonzero_count / edge_count,
      passesTolerance = if (edge_count == 0) TRUE else exceedance_count / edge_count <= 0.05
    ),
    tolerance = tolerance,
    firstLabel = first_label,
    secondLabel = second_label
  )
}

#' Compare directed bootstrap strengths across backends
#'
#' Quantifies how similar two learned structures are at the directed-edge bootstrap
#' level. Missing edges are treated as zero strength in the other structure.
#'
#' @param firstStructure A result from `LearnBayesianNetwork()`.
#' @param secondStructure A result from `LearnBayesianNetwork()`.
#' @param firstLabel Label for the first structure.
#' @param secondLabel Label for the second structure.
#' @param tolerance Absolute difference threshold used to count exceeded edges.
#'
#' @return A list with `$edgeTable`, `$summary`, `$tolerance`, `$firstLabel`, and `$secondLabel`.
CompareBackendBootstrapStrengths <- function(firstStructure, secondStructure, firstLabel = "First", secondLabel = "Second", tolerance = 0.05) {
  if (!is.numeric(tolerance) || length(tolerance) != 1 || is.na(tolerance) || tolerance < 0) {
    stop("tolerance must be a single non-negative numeric value.", call. = FALSE)
  }

  compare_directed_bootstrap_strength(
    firstStructure$bootstrapStrength,
    secondStructure$bootstrapStrength,
    firstLabel,
    secondLabel,
    tolerance
  )
}

# Build a per-edge stage-probability comparison table from two sevt models.
# Returns a data frame with columns: edge, firstTV, secondTV, deltaTV, direction.
build_stage_comparison_table <- function(first_sevt, second_sevt, first_label, second_label) {
  assoc_first  <- compute_edge_association(first_sevt)
  assoc_second <- compute_edge_association(second_sevt)

  all_edges <- union(names(assoc_first), names(assoc_second))
  if (length(all_edges) == 0) {
    return(data.frame(
      edge = character(0), firstStrength = numeric(0), secondStrength = numeric(0),
      firstSignedStrength = numeric(0), secondSignedStrength = numeric(0),
      delta = numeric(0), absoluteDelta = numeric(0), direction = character(0),
      stringsAsFactors = FALSE
    ))
  }

  first_tv  <- vapply(all_edges, function(e) { v <- assoc_first[e];  if (length(v) == 0 || is.na(v)) 0 else v }, numeric(1))
  second_tv <- vapply(all_edges, function(e) { v <- assoc_second[e]; if (length(v) == 0 || is.na(v)) 0 else v }, numeric(1))
  delta     <- first_tv - second_tv

  data.frame(
    edge                = all_edges,
    firstStrength       = first_tv,
    secondStrength      = second_tv,
    firstSignedStrength = first_tv,
    secondSignedStrength = second_tv,
    delta               = delta,
    absoluteDelta       = abs(delta),
    direction           = dplyr::case_when(
                            delta >  0.1 ~ paste0("Stronger in ", first_label),
                            delta < -0.1 ~ paste0("Stronger in ", second_label),
                            TRUE         ~ "Similar"
                          ),
    stringsAsFactors = FALSE
  ) |>
    dplyr::arrange(dplyr::desc(.data$absoluteDelta))
}

#' Compare two staged event tree structures
#'
#' Computes stage-probability-based association differences between two stagedtrees
#' structure results, along with Hamming distance and CID metrics.
#'
#' @param firstStructure A result from `LearnBayesianNetwork()` with a stagedtrees backend.
#' @param secondStructure A result from `LearnBayesianNetwork()` with a stagedtrees backend.
#' @param firstLabel Label for the first structure.
#' @param secondLabel Label for the second structure.
#'
#' @return A list with `$comparisonTable`, `$firstLabel`, `$secondLabel`,
#'   `$hammingDistance`, and `$cid`.
CompareStages <- function(firstStructure, secondStructure, firstLabel = "First", secondLabel = "Second") {
  if (is.null(firstStructure$sevtModel) || is.null(secondStructure$sevtModel)) {
    stop("CompareStages() requires both structures to have a $sevtModel (stagedtrees backend).", call. = FALSE)
  }

  first_sevt  <- firstStructure$sevtModel
  second_sevt <- secondStructure$sevtModel

  comparison_table <- build_stage_comparison_table(first_sevt, second_sevt, firstLabel, secondLabel)

  hamming_dist <- NA_integer_
  cid_value    <- NA_real_

  if (requireNamespace("clue", quietly = TRUE)) {
    hamming_dist <- tryCatch(
      stagedtrees::hamming_stages(first_sevt, second_sevt),
      error = function(e) NA_integer_
    )
  } else {
    message("Install the 'clue' package to enable Hamming distance comparison.")
  }

  cid_result <- tryCatch(
    stagedtrees::cid(first_sevt, second_sevt),
    error = function(e) list(cid = NA_real_)
  )
  cid_value <- cid_result$cid %||% NA_real_

  list(
    comparisonTable = comparison_table,
    firstLabel      = firstLabel,
    secondLabel     = secondLabel,
    hammingDistance = hamming_dist,
    cid             = cid_value
  )
}

#' Compare two learned network structures
#'
#' Dispatches to `CompareStages()` when both results use the stagedtrees backend,
#' otherwise computes signed bootstrap-strength differences (bnlearn backend).
#'
#' @param firstStructure A result from `LearnBayesianNetwork()`.
#' @param secondStructure A result from `LearnBayesianNetwork()`.
#' @param firstLabel Label for the first structure.
#' @param secondLabel Label for the second structure.
#'
#' @return A list containing a comparison table and summary metrics.
CompareStructures <- function(firstStructure, secondStructure, firstLabel = "First", secondLabel = "Second") {
  both_stagedtrees <- !is.null(firstStructure$sevtModel) && !is.null(secondStructure$sevtModel)

  if (both_stagedtrees) {
    return(CompareStages(firstStructure, secondStructure, firstLabel, secondLabel))
  }

  # bnlearn signed bootstrap strength path
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
      delta         = .data$firstSignedStrength - .data$secondSignedStrength,
      absoluteDelta = abs(.data$delta),
      direction     = dplyr::case_when(
                        .data$delta >  0.1 ~ paste0("Stronger in ", firstLabel),
                        .data$delta < -0.1 ~ paste0("Stronger in ", secondLabel),
                        TRUE               ~ "Similar"
                      ),
      edgeType = vapply(.data$edge, function(edge_name) {
        edge_nodes <- strsplit(edge_name, " -- ", fixed = TRUE)[[1]]
        classify_edge_type(edge_nodes[[1]], edge_nodes[[2]], active_exogenous_variables, measurement_groups)
      }, character(1))
    ) |>
    dplyr::arrange(dplyr::desc(.data$absoluteDelta))

  list(
    comparisonTable = comparison_table,
    firstLabel      = firstLabel,
    secondLabel     = secondLabel,
    hammingDistance = NA_integer_,
    cid             = NA_real_
  )
}
