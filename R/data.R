#' Simple simulated two-sex model (r4ss "simple_small"): SS3 report
#'
#' Output of \code{r4ss::SS_output()} for the "simple_small" Stock Synthesis (3.30.26) example of the r4ss package,
#' the example also used by the PIFSC ss3diags package (\url{https://github.com/PIFSCstockassessments/ss3diags}), so
#' results can be compared between ss3diags and ss3diags2. Simulated stock, 2011-2022, one season, one area, two
#' sexes, one fishery (FISHERY: catch, length and age compositions) and two surveys (SURVEY1, SURVEY2: indices, length
#' and age compositions); Hessian available (CoVar reduced to the derived quantities), initial F, Beverton-Holt
#' stock-recruitment. The SS3 input and output files are in \code{system.file("extdata", "simple_small", package = "r4ss")}.
#'
#' Further examples: the two cookbook models (\code{\link{ss3sma}}, \code{\link{ss3phk}}) with retrospective runs.
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

#' Cookbook example: North Atlantic shortfin mako (Isurus oxyrinchus)
#'
#' SS3 reports of the North Atlantic shortfin mako assessment used in the diagnostics cookbook (Carvalho et al.
#' 2021) and in ss3diags: \code{ss3sma} the reference model, \code{retro.sma} the retrospective runs (list of five
#' peels plus the reference, as \code{r4ss::SSgetoutput()}), \code{aspm.sma} the age-structured production model
#' (summary for model comparison). Copied from the ss3diags package (\code{natl.sma}); CoVar reduced to the derived
#' quantities and large age/length tables dropped.
#'
#' @format lists (r4ss output)
#' @source \url{https://github.com/jabbamodel/ss3diags}
#' @references Carvalho, F., Winker, H., Courtney, D., Kapur, M., Kell, L., Cardinale, M., Schirripa, M.,
#'   Kitakado, T., Yemane, D., Piner, K.R., Maunder, M.N., Taylor, I., Wetzel, C.R., Doering, K., Johnson, K.F. and
#'   Methot, R.D. (2021). A cookbook for using model diagnostics in integrated stock assessments. Fisheries
#'   Research 240: 105959.
#' @examples
#' data(ss3sma); data(retro.sma)
#' ss3_runs(ss3sma)
#' ss3_mohn(ss3_retro(retro.sma), c("SSB", "F/Ftgt"))
#' @name sma
#' @aliases ss3sma retro.sma aspm.sma
NULL

#' @rdname sma
"ss3sma"
#' @rdname sma
"retro.sma"
#' @rdname sma
"aspm.sma"

#' Cookbook example: Pacific hake (Merluccius productus)
#'
#' SS3 reports of the Pacific hake assessment used in the diagnostics cookbook (Carvalho et al. 2021) and in
#' ss3diags: \code{ss3phk} the reference model (Dirichlet-multinomial age compositions, empirical weight at age),
#' \code{retro.phk} the retrospective runs, \code{aspm.phk} the age-structured production model. Copied from the
#' ss3diags package (\code{pac.hke}), reduced as for \code{\link{ss3sma}}.
#'
#' @format lists (r4ss output)
#' @source \url{https://github.com/jabbamodel/ss3diags}
#' @references Carvalho et al. (2021), see \code{\link{ss3sma}}.
#' @examples
#' data(ss3phk); data(retro.phk)
#' hc <- ss3_hindcast(retro.phk)
#' hc$mase
#' @name phk
#' @aliases ss3phk retro.phk aspm.phk
NULL

#' @rdname phk
"ss3phk"
#' @rdname phk
"retro.phk"
#' @rdname phk
"aspm.phk"
