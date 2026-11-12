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
#' equilibrium yield table is empty: \code{\link{ss3stk}()} then builds the trajectories from the time series
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

#' Simple simulated two-sex model (r4ss "simple_small"): SS3 report
#'
#' Output of \code{r4ss::SS_output()} for the "simple_small" Stock Synthesis (3.30.26) example of the r4ss package,
#' the example also used by the PIFSC ss3diags package (\url{https://github.com/PIFSCstockassessments/ss3diags}), so
#' results can be compared between ss3diags and ss3diags2. Simulated stock, 2011-2022, one season, one area, two
#' sexes, one fishery (FISHERY: catch, length and age compositions) and two surveys (SURVEY1, SURVEY2: indices, length
#' and age compositions); Hessian available (CoVar reduced to the derived quantities), initial F, Beverton-Holt
#' stock-recruitment. The SS3 input and output files are in \code{system.file("extdata", "simple_small", package = "r4ss")}.
#'
#' A season-, sex- and area-structured example (northern shrimp pra.27.3a4a: four seasons, hermaphrodite, two areas,
#' six area-specific fleets) is shipped as an rds file with its SS3 input files:
#' \code{readRDS(system.file("extdata", "pra.27.3a4a", "pra.27.3a4a.rds", package = "ss3diags2"))}; see the
#' vignette "Season-, sex- and area-structured models".
#'
#' @format list (r4ss \code{SS_output()})
#' @source r4ss package (\url{https://github.com/r4ss/r4ss}), \code{inst/extdata/simple_small}
#' @examples
#' data(ss3simple)
#' ss3_fleets(ss3simple)
#' stk <- ss3stk(ss3simple)
"ss3simple"
