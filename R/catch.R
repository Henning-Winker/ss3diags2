# catch.R - catch by fleet

.catch_cols <- c(obs = "Obs", exp = "Exp", ret_bio = "ret_bio", kill_bio = "kill_bio", sel_bio = "sel_bio",
                 ret_num = "ret_num", kill_num = "kill_num", sel_num = "sel_num", F = "F")

#' Catch by fleet
#'
#' Extracts the catch of each SS3 catch fleet (fleet types 1 and 2) from \code{ss3rep$catch} as \code{FLQuants}, one
#' \code{FLQuant} per fleet (quant \code{"all"}, season and area dimensions as in the model), or as a long table.
#' This is the catch input of JABBA3 (\code{FLQuants} by fleet) and the relative catch by gear of FLicc
#' (\code{catch_by_gear}).
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param fleets fleet numbers or names (default: all catch fleets with catch)
#' @param what \code{"obs"} (observed catch, default), \code{"exp"} (expected), model-estimated catch by type
#'   (\code{"ret_bio"}, \code{"kill_bio"}, \code{"sel_bio"}: biomass retained, killed, selected;
#'   \code{"ret_num"}, \code{"kill_num"}, \code{"sel_num"}: the same in numbers) or the fleet \code{"F"}
#' @param years optional years to keep (default: start to end year of the model)
#' @param long return the long data.table instead of FLQuants
#' @return FLQuants (named by fleet) or a data.table with columns \code{fleet}, \code{name}, \code{age} ("all"),
#'   \code{year}, \code{unit}, \code{season}, \code{area}, \code{iter}, \code{data}, \code{obs}, \code{exp},
#'   \code{se} (log-scale catch SE) and \code{units}. Units: "t" for biomass, "1000" for numbers (SS3 convention),
#'   "f" for F.
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' ss3_catch(ss3sma)                       # observed catch, FLQuants by fleet
#' ss3_catch(ss3sma, what = "kill_bio")    # catch in weight from the model, for fleets reported in numbers
#' }
#' @export
ss3_catch <- function(ss3rep, fleets = NULL, what = c("obs", "exp", "ret_bio", "kill_bio", "sel_bio", "ret_num",
                      "kill_num", "sel_num", "F"), years = NULL, long = FALSE){
  what <- match.arg(what)
  fl <- ss3_fleets(ss3rep)
  ct <- as.data.table(ss3rep$catch)
  if(!nrow(ct)) stop("ss3rep$catch is empty", call. = FALSE)
  yrs <- years %||% seq(ss3rep$startyr, ss3rep$endyr)
  ct <- ct[Yr %in% yrs]
  avail <- fl[type %in% c("catch", "bycatch") & fleet %in% unique(ct$Fleet)]$name
  fleets <- .pick_fleets(ss3rep, fleets, avail)
  ns <- ss3rep$nseasons %||% 1; na <- ss3rep$nareas %||% 1
  col <- .catch_cols[[what]]
  if(!col %in% names(ct)) stop("column '", col, "' not found in ss3rep$catch", call. = FALSE)
  d <- data.table(fleet = ct$Fleet, name = .fleetnames(ss3rep)[ct$Fleet], age = "all", year = ct$Yr,
                  unit = "unique", season = .season_lab(ct$Seas, ns), area = .area_lab(ct$Area, na), iter = "1",
                  data = as.numeric(ct[[col]]), obs = ct$Obs, exp = ct$Exp, se = .col(ct, c("se", "SE"), NA_real_))
  d <- d[name %in% fleets]
  cu <- fl$catch_units[match(d$fleet, fl$fleet)]
  d[, units := switch(what, obs = , exp = ifelse(cu == "numbers", "1000", "t"), F = "f",
                      if(grepl("_num$", what)) "1000" else "t")]
  # several catch records within the same cell (e.g. months within a season) are summed
  if(anyDuplicated(d, by = c("name", "year", "season", "area")))
    d <- d[, .(fleet = fleet[1], age = "all", unit = "unique", iter = "1", data = sum(data), obs = sum(obs),
               exp = sum(exp), se = mean(se), units = units[1]), by = .(name, year, season, area)]
  setcolorder(d, c("fleet", "name", "age", "year", "unit", "season", "area", "iter", "data"))
  if(long) return(d[])
  u <- d[, .(u = units[1]), by = name]
  long2flqs(d, quant = "age", units = u$u[match(fleets, u$name)], order = fleets, years = yrs)
}
