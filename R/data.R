#' Synthetic three-level landscape
#'
#' A 20 x 15 grid of 1 km cells (EPSG:3035) with planted structure:
#' three regions (signal in `temp`), two local types shared across regions
#' (signal in `prec`), and region-specific local-type proportions.
#'
#' @format An `sf` data frame with 300 rows:
#' \describe{
#'   \item{ID}{Integer cell ID.}
#'   \item{temp}{Simulated temperature; carries the region signal.}
#'   \item{prec}{Simulated precipitation; carries the local-type signal.}
#'   \item{true_region, true_local}{Planted labels. Exclude from analysis variables.}
#' }
#' @source Simulated; see `data-raw/pulse_grid.R`.
"pulse_grid"