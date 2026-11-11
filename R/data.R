#' Black anglerfish (Lophius budegassa) in ICES areas 7 and 8abd: SS3 report
#'
#' Output of \code{r4ss::SS_output()} for the Stock Synthesis (3.30.22) assessment of black anglerfish in ICES
#' divisions 7 and 8.a, b, d (ank.27.78abd), 1950-2025. One season (four sub-seasons), two sexes, one area, four fleets:
#' two fishing fleets (FL1_Trawls, FL2_Gillnets) with retention (discard and retained length compositions, observed
#' discards) and two surveys (FR_IE_IBTS, IE_IAMS) with sex-specific length compositions and biomass indices.
#' The SS3 input files (data.ss, control.ss, starter.ss, forecast.ss) are in
#' \code{system.file("extdata", "ank.27.78abd", package = "ss3diags2")}.
#'
#' The run was read without CoVar and its derived quantities are reported as 0 (no hessian / sdreport), and the
#' equilibrium yield table is empty: \code{\link{ss3_FLStockR}()} then builds the trajectories from the time series
#' (\code{source = "timeseries"}), MVLN draws, Kobe plots and production functions are not available.
#'
#' @format list (r4ss \code{SS_output()})
#' @source ICES WGBIE assessment of black anglerfish (ank.27.78abd)
#' @examples
#' data(ss3ank)
#' ss3_fleets(ss3ank)
#' lfd <- ss3_lencomp(ss3ank)
#' names(lfd)
"ss3ank"
