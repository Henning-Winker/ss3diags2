# fleets.R - SS3 fleet table

.fleet_types  <- c("1" = "catch", "2" = "bycatch", "3" = "survey", "4" = "predator")
.catch_units  <- c("1" = "biomass", "2" = "numbers")
.survey_units <- c("0" = "numbers", "1" = "biomass", "2" = "F", "30" = "spawnbio", "31" = "recdev", "32" = "spawnbio",
                   "33" = "recruitment", "34" = "depletion", "35" = "parmdev", "36" = "recdev")

#' SS3 fleet table
#'
#' One row per SS3 fleet with its type, area, units and the number of observations of each data type, used to
#' select fleets for the extraction functions and to map them to JABBA3 (catch fleets, index fleets) and FLicc
#' (gears with length compositions).
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @return data.table with columns \code{fleet}, \code{name}, \code{type} (catch, bycatch, survey, predator),
#'   \code{area}, \code{catch_units} (biomass or numbers; catch in numbers is in 1000s), \code{survey_units}
#'   (numbers, biomass, F, ...), and counts \code{ncatch} (years with catch > 0), \code{nindex}, \code{nlen},
#'   \code{nage} (number of composition samples)
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' ss3_fleets(ss3sma)
#' }
#' @export
ss3_fleets <- function(ss3rep){
  nm <- .fleetnames(ss3rep); nf <- length(nm); id <- seq_len(nf)
  def <- ss3rep$definitions
  pick <- function(vec, col, default){
    v <- ss3rep[[vec]]
    if(is.null(v) && !is.null(def) && col %in% names(def)) v <- def[[col]]
    if(is.null(v)) rep(default, nf) else rep_len(v, nf)
  }
  ft <- pick("fleet_type", "fleet_type", 1)
  ar <- pick("fleet_area", "area", 1)
  cu <- pick("catch_units", "catch_units", 1)
  su <- pick("survey_units", "survey_units", 1)
  cnt <- function(db, f){
    if(is.null(db) || !NROW(db)) return(0L)
    d <- as.data.table(db)[Fleet == f]
    if(!nrow(d)) return(0L)
    nrow(unique(d[, intersect(c("Yr", "Seas", "Month", "Sexes", "Part", "Ageerr"), names(d)), with = FALSE]))
  }
  ct <- as.data.table(ss3rep$catch)
  data.table(fleet = id, name = nm,
             type = unname(.fleet_types[as.character(ft)]),
             area = ar,
             catch_units = ifelse(ft %in% 1:2, unname(.catch_units[as.character(cu)]), NA_character_),
             survey_units = unname(.survey_units[as.character(su)]),
             ncatch = vapply(id, function(f) if(nrow(ct)) sum(ct$Fleet == f & ct$Obs > 0, na.rm = TRUE) else 0L, 1L),
             nindex = vapply(id, function(f) if(NROW(ss3rep$cpue)) sum(ss3rep$cpue$Fleet == f) else 0L, 1L),
             nlen = vapply(id, function(f) cnt(ss3rep$lendbase, f), 1L),
             nage = vapply(id, function(f) cnt(ss3rep$agedbase, f), 1L))
}

# resolve fleets given by number or name to names, restricted to 'avail'
.pick_fleets <- function(ss3rep, fleets, avail){
  nm <- .fleetnames(ss3rep)
  if(is.null(fleets)) return(intersect(nm, avail))
  if(is.numeric(fleets)) fleets <- nm[fleets]
  bad <- setdiff(fleets, avail)
  if(length(bad)) warning("no data for fleet(s): ", paste(bad, collapse = ", "), call. = FALSE)
  intersect(fleets, avail)
}
