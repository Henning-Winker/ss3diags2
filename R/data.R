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

#' Northern shrimp (Pandalus borealis) in Skagerrak and the Norwegian Deep: SS3 report
#'
#' Output of \code{r4ss::SS_output()} for the Stock Synthesis (3.30.22) assessment of northern shrimp in ICES
#' division 3.a and subarea 4.a East (pra.27.3a4a), 1970-2024. A season-, sex- and area-structured model:
#' \itemize{
#'   \item four seasons (quarters), spawning in season 1;
#'   \item two sexes with protandric hermaphroditism (born male, change to female; unit "M" and "F");
#'   \item two areas (1: Skagerrak 3.a, 2: Norwegian Deep 4.a) with time-varying recruitment distribution;
#'   \item six area-specific fleets (DK3a, DK4a, NOR3a, NOR4a, SWE3a, SWE4a) with seasonal catches, retention and
#'     discard, retained and combined length compositions, and two surveys (NORSURVEY3a, NORSURVEY4a) with
#'     sex-specific length compositions and biomass indices;
#'   \item F reported as Fbar over ages 1-3 (Z - M), MSY proxy F at B30%.
#' }
#' The CoVar table is reduced to the derived quantities and the length-based numbers at length are dropped to keep
#' the package small. The SS3 input files are in \code{system.file("extdata", "pra.27.3a4a", package = "ss3diags2")}.
#'
#' @format list (r4ss \code{SS_output()})
#' @source ICES assessment of northern shrimp (pra.27.3a4a)
#' @examples
#' data(ss3pra)
#' ss3_fleets(ss3pra)
#' pf <- ss3_partialF(ss3pra, groups = "area", method = "fbar")
#' imp <- ss3_impact(ss3pra, groups = list(Denmark = 1:2, Norway = 3:4, Sweden = 5:6))
"ss3pra"
