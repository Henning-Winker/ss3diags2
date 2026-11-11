# models.R - SS3 observations shaped as inputs of JABBA3 and FLicc

# catch in weight for each fleet: observed where reported in weight, model retained/killed weight otherwise
.catch_weight <- function(ss3rep, fleets, catch, years){
  fl <- ss3_fleets(ss3rep)
  if(catch != "auto") return(ss3_catch(ss3rep, fleets, what = catch, years = years, long = TRUE))
  d <- ss3_catch(ss3rep, fleets, what = "obs", years = years, long = TRUE)
  num <- unique(d$name[fl$catch_units[match(d$fleet, fl$fleet)] == "numbers"])
  if(length(num)){
    message("catch in numbers for ", paste(num, collapse = ", "), ": using the SS3 retained catch in weight (ret_bio)")
    w <- ss3_catch(ss3rep, num, what = "ret_bio", years = years, long = TRUE)
    d <- rbind(d[!name %in% num], w)
  }
  d[order(match(name, .fleetnames(ss3rep)), year, season)]
}

# sum seasons (catch) or average them (indices, ignoring NA) to an annual FLQuant
.annual <- function(x, fun = c("sum", "mean")){
  fun <- match.arg(fun)
  if(dim(x)[4] == 1) return(x)
  if(fun == "sum") FLCore::seasonSums(x, na.rm = TRUE) else FLCore::seasonMeans(x, na.rm = TRUE)
}

#' SS3 observations as JABBA3 inputs
#'
#' Extracts catch by fleet and abundance indices from an SS3 report in the form \code{JABBA3::jb3_data()} takes them:
#' \code{catch} as \code{FLQuants} by fleet (annual, in weight) and \code{indices} as \code{FLIndices} of
#' \code{FLIndexBiomass} (annual, input log-scale SE in \code{index.var}, as JABBA3 estimates an additive SD on top).
#'
#' JABBA3 is a biomass dynamic model, so catch must be in weight. With \code{catch = "auto"} (default) fleets with
#' catch in weight keep the observed catch and fleets with catch in numbers get the SS3 estimate of retained catch in
#' weight (\code{ret_bio}), with a message. Seasons are summed (catch) or averaged (indices) for annual models.
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param catch.fleets,index.fleets fleets to keep (numbers or names; default all catch fleets and all index fleets)
#' @param catch \code{"auto"} (default), \code{"obs"}, \code{"ret_bio"} or \code{"kill_bio"}
#' @param se \code{"input"} (default) or \code{"total"} index SE, see \code{\link{ss3_indices}()}
#' @param annual collapse seasons to years (default TRUE; JABBA3 seasonal models take seasonal data frames)
#' @param combine.catch also return the total catch over fleets as a single FLQuant (\code{catch.total})
#' @return list with \code{catch} (FLQuants), \code{indices} (FLIndices), \code{catch.total} (FLQuant) and
#'   \code{fleets} (the \code{\link{ss3_fleets}()} table)
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' inp <- ss2jabba3(ss3sma)
#' # d <- JABBA3::jb3_data(inp$catch, inp$indices, stock = "sma")
#' }
#' @export
ss2jabba3 <- function(ss3rep, catch.fleets = NULL, index.fleets = NULL, catch = c("auto", "obs", "ret_bio", "kill_bio"),
                      se = c("input", "total"), annual = TRUE, combine.catch = TRUE){
  catch <- match.arg(catch); se <- match.arg(se)
  yrs <- seq(ss3rep$startyr, ss3rep$endyr)
  d <- .catch_weight(ss3rep, catch.fleets, catch, yrs)
  d <- d[, .(fleet = fleet[1], data = sum(data), units = units[1]), by = .(name, age, year, unit, season, iter,
                                                                            area)]
  if(annual) d <- d[, .(fleet = fleet[1], season = "all", data = sum(data, na.rm = TRUE), units = units[1]),
                    by = .(name, age, year, unit, iter, area)]
  # JABBA3 has no area dimension: sum areas
  d <- d[, .(fleet = fleet[1], area = "unique", data = sum(data, na.rm = TRUE), units = units[1]),
         by = .(name, age, year, unit, season, iter)]
  nm <- unique(d$name)
  ct <- long2flqs(d, quant = "age", units = "t", order = nm, years = yrs)
  idx <- ss3_indices(ss3rep, fleets = index.fleets, se = se)
  if(annual) idx <- FLCore::FLIndices(lapply(idx, function(x){
    if(dim(FLCore::index(x))[4] == 1) return(x)
    y <- FLCore::FLIndexBiomass(index = .annual(FLCore::index(x), "mean"), name = FLCore::name(x),
                                desc = FLCore::desc(x))
    FLCore::index.var(y) <- .annual(FLCore::index.var(x), "mean")
    FLCore::index.q(y) <- .annual(FLCore::index.q(x), "mean")
    y@range[c("startf", "endf")] <- x@range[c("startf", "endf")]
    y@distribution <- x@distribution
    y
  }))
  out <- list(catch = ct, indices = idx, fleets = ss3_fleets(ss3rep))
  if(combine.catch) out$catch.total <- Reduce(`+`, lapply(ct, function(q){ q[is.na(q)] <- 0; q }))
  out
}

#' SS3 observations as FLicc inputs
#'
#' Extracts length compositions by fleet (gear) and the matching catch by gear from an SS3 report in the form
#' \code{FLicc} takes them: \code{lfd} as \code{FLQuants}, one per gear, of numbers at length (dimensions
#' \code{len} x \code{year}, length bins as lower bounds in cm, sexes and seasons summed, as \code{FLicc::FLQuantLen()}
#' with \code{midL = FALSE}), all gears on the same length bins and years, and \code{catch_by_gear} as
#' \code{FLQuants} of catch in weight for the same gears.
#'
#' FLicc fits all gears in each year, so with \code{years = "common"} only years in which every gear has a length
#' sample are kept; with \code{years = "all"} gears without a sample in a year are \code{NA}.
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param fleets gears to keep (numbers or names; default: catch fleets with length compositions)
#' @param years \code{"common"} (default), \code{"all"} or a numeric vector of years
#' @param nsamp sample size used to convert SS3 proportions to numbers, \code{"adj"} (default) or \code{"input"}
#' @param catch catch in weight by gear, as in \code{\link{ss2jabba3}()}
#' @param surveys also include survey fleets with length compositions (default FALSE: FLicc gears are catch fleets)
#' @return list with \code{lfd} (FLQuants), \code{catch_by_gear} (FLQuants, annual, in weight; NULL for surveys)
#'   and \code{years}
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' inp <- ss2flicc(ss3sma, years = "all")
#' # stk <- FLicc::stocklen(inp$lfd, lhpar)
#' }
#' @export
ss2flicc <- function(ss3rep, fleets = NULL, years = "common", nsamp = c("adj", "input"),
                     catch = c("auto", "obs", "ret_bio", "kill_bio"), surveys = FALSE){
  nsamp <- match.arg(nsamp); catch <- match.arg(catch)
  fl <- ss3_fleets(ss3rep)
  avail <- fl[nlen > 0 & (surveys | type %in% c("catch", "bycatch"))]$name
  fleets <- .pick_fleets(ss3rep, fleets, avail)
  if(!length(fleets)) stop("no fleets with length compositions", call. = FALSE)
  d <- ss3_lencomp(ss3rep, fleets = fleets, scale = "n", nsamp = nsamp, sex = "combine", long = TRUE)
  # FLicc: one composition per gear and year (sum partitions, seasons and areas)
  d <- d[, .(data = sum(data, na.rm = TRUE)), by = .(name, len, year)]
  yobs <- d[, .(yrs = list(sort(unique(year)))), by = name]
  yrs <- if(is.numeric(years)) years else if(years == "common") Reduce(intersect, yobs$yrs) else
    sort(unique(unlist(yobs$yrs)))
  if(!length(yrs)) stop("no year with length samples for all gears; use years = \"all\"", call. = FALSE)
  # common length bins: all data bins of the model
  lens <- sort(unique(d$len))
  grid <- CJ(name = fleets, len = lens, year = yrs)
  d <- d[grid, on = .(name, len, year)]
  d[is.na(data) & paste(name, year) %in% d[!is.na(data), paste(name, year)], data := 0]
  lfd <- FLCore::FLQuants(stats::setNames(lapply(fleets, function(f){
    q <- long2flq(d[name == f], quant = "len", units = "cm", years = yrs)
    q[, as.character(yrs)]
  }), fleets))
  cg <- intersect(fleets, fl[type %in% c("catch", "bycatch")]$name)
  cbg <- NULL
  if(length(cg)){
    cw <- .catch_weight(ss3rep, cg, catch, seq(ss3rep$startyr, ss3rep$endyr))
    cw <- cw[, .(age = "all", data = sum(data, na.rm = TRUE)), by = .(name, year)]
    cbg <- long2flqs(cw[year %in% yrs], quant = "age", units = "t", order = cg, years = yrs)
  }
  list(lfd = lfd, catch_by_gear = cbg, years = yrs)
}
