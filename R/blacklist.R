build_variable_order <- function(active_exogenous_variables, active_measurement_groups) {
  # Sort measurement groups by tier ascending, then alphabetically within each tier
  tier_values <- vapply(active_measurement_groups, `[[`, numeric(1), "tier")
  sorted_groups <- active_measurement_groups[order(tier_values)]

  measurement_columns_ordered <- unlist(lapply(sorted_groups, function(group) {
    sort(group$columns)
  }))

  c(sort(active_exogenous_variables), measurement_columns_ordered)
}

build_measurement_to_exogenous_blacklist <- function(measurement_columns, active_exogenous_variables) {
  if (length(measurement_columns) == 0 || length(active_exogenous_variables) == 0) {
    return(data.frame(from = character(0), to = character(0)))
  }

  expand.grid(from = measurement_columns, to = active_exogenous_variables, stringsAsFactors = FALSE)
}

build_exogenous_to_exogenous_blacklist <- function(active_exogenous_variables) {
  if (length(active_exogenous_variables) < 2) {
    return(data.frame(from = character(0), to = character(0)))
  }

  expand.grid(from = active_exogenous_variables, to = active_exogenous_variables, stringsAsFactors = FALSE) |>
    dplyr::filter(.data$from != .data$to)
}

build_tier_blacklist <- function(measurement_groups) {
  blacklist_rows <- list()
  row_index <- 1

  for (source_group in measurement_groups) {
    for (target_group in measurement_groups) {
      if (source_group$tier > target_group$tier && length(source_group$columns) > 0 && length(target_group$columns) > 0) {
        blacklist_rows[[row_index]] <- expand.grid(
          from = source_group$columns,
          to = target_group$columns,
          stringsAsFactors = FALSE
        )
        row_index <- row_index + 1
      }
    }
  }

  if (length(blacklist_rows) == 0) {
    return(data.frame(from = character(0), to = character(0)))
  }

  dplyr::bind_rows(blacklist_rows)
}

#' Build a blacklist from variable tiers and study design
#'
#' @param discretizationResult A result from `DiscretizeMeasurements()`.
#' @param configuration A networkR configuration object.
#'
#' @return A list containing the blacklist and active exogenous variables.
BuildBlacklist <- function(discretizationResult, configuration) {
  discrete_data <- discretizationResult$discreteData
  active_exogenous_variables <- find_active_exogenous_variables(
    discrete_data,
    discretizationResult$exogenousVariables
  )
  active_measurement_columns <- find_active_measurement_columns(
    discrete_data,
    discretizationResult$measurementColumns
  )
  active_measurement_groups <- filter_measurement_groups(
    discretizationResult$measurementGroups,
    active_measurement_columns
  )

  blacklist <- dplyr::bind_rows(
    build_measurement_to_exogenous_blacklist(
      active_measurement_columns,
      active_exogenous_variables
    ),
    build_exogenous_to_exogenous_blacklist(active_exogenous_variables),
    build_tier_blacklist(active_measurement_groups)
  ) |>
    dplyr::distinct(.data$from, .data$to)

  list(
    blacklist = blacklist,
    variableOrder = build_variable_order(active_exogenous_variables, active_measurement_groups),
    activeExogenousVariables = active_exogenous_variables,
    measurementGroups = active_measurement_groups,
    measurementColumns = active_measurement_columns
  )
}