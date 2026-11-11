# utils.R - internal helpers: SS3 dimensions, labels and long table -> FLQuant conversion

`%||%` <- function(a, b) if(is.null(a)) b else a

.need <- function(pkg, msg = paste0("package '", pkg, "' is required")){
  if(!requireNamespace(pkg, quietly = TRUE)) stop(msg, call. = FALSE)
  invisible(TRUE)
}

# first column of x that exists, else NULL (r4ss renames columns across versions)
.col <- function(x, cand, default = NULL){
  hit <- cand[cand %in% names(x)]
  if(length(hit)) x[[hit[1]]] else default
}

# FLQuant dimension columns, in order
.flqdims <- c("year", "unit", "season", "area", "iter")

.season_lab <- function(seas, ns) if(ns == 1) rep("all", length(seas)) else as.character(seas)
.area_lab   <- function(area, na) if(na == 1) rep("unique", length(area)) else as.character(area)

# sex as FLR unit: "unique" for combined-sex data (SS3 data option Sexes = 0) or single-sex models, else F/M
.unit_lab <- function(sex, sexes, nsexes){
  if(nsexes == 1) return(rep("unique", length(sex)))
  ifelse(sexes == 0, "unique", c("F", "M")[sex])
}

# SS3 fleet names in model order
.fleetnames <- function(ss3rep) as.character(ss3rep$FleetNames)

# order labels: numeric if possible, FLR keywords first
.lev <- function(x, first = c("unique", "all", "F", "M")){
  u <- unique(as.character(x))
  num <- suppressWarnings(as.numeric(u))
  if(all(!is.na(num))) return(u[order(num)])
  c(intersect(first, u), sort(setdiff(u, first)))
}

#' Convert a long table to an FLQuant
#'
#' Builds an \code{FLQuant} from a long table with columns \code{quant} (e.g. \code{len} or \code{age}), \code{year},
#' \code{unit}, \code{season}, \code{area}, \code{iter} and \code{data}, keeping numeric order of length bins and
#' ages and filling years without data with \code{NA} (FLR convention: contiguous years). Missing dimension columns
#' are filled with the FLR defaults (\code{"all"}, \code{"unique"}, 1).
#'
#' @param x data.frame or data.table
#' @param quant name of the first dimension column; \code{"age"} with value \code{"all"} for aggregated quantities
#' @param units units of the FLQuant
#' @param years optional year range (default: min to max of \code{x$year})
#' @param fill value for cells without data (default NA)
#' @return FLQuant
#' @examples
#' d <- data.frame(len = c(5, 10, 15, 5, 10, 15), year = rep(c(2000, 2002), each = 3), data = 1:6)
#' long2flq(d, quant = "len", units = "numbers")
#' @export
long2flq <- function(x, quant = "age", units = "NA", years = NULL, fill = NA_real_){
  x <- as.data.table(x)
  if(!quant %in% names(x)) x[, (quant) := "all"]
  defs <- list(unit = "unique", season = "all", area = "unique", iter = "1")
  for(d in names(defs)) if(!d %in% names(x)) x[, (d) := defs[[d]]]
  if(!"year" %in% names(x)) stop("long2flq: column 'year' required", call. = FALSE)
  dims <- c(quant, .flqdims)
  if(anyDuplicated(x, by = dims)) stop("long2flq: duplicated rows for ", paste(dims, collapse = ", "),
                                       "; aggregate first", call. = FALSE)
  yrs <- years %||% seq(min(x$year), max(x$year))
  dmn <- list(.lev(x[[quant]]), as.character(yrs), .lev(x$unit), .lev(x$season), .lev(x$area), .lev(x$iter))
  names(dmn) <- dims
  arr <- array(fill, dim = lengths(dmn), dimnames = dmn)
  x <- x[as.character(year) %in% dmn$year]
  idx <- vapply(seq_along(dims), function(i) match(as.character(x[[dims[i]]]), dmn[[i]]), integer(nrow(x)))
  arr[matrix(idx, ncol = length(dims))] <- x$data
  FLCore::FLQuant(arr, units = units)
}

#' @rdname long2flq
#' @param by column defining the FLQuants elements (default \code{"name"}, the fleet name)
#' @param order optional order of the elements
#' @export
long2flqs <- function(x, quant = "age", units = "NA", by = "name", order = NULL, years = NULL){
  x <- as.data.table(x)
  nm <- order %||% unique(x[[by]])
  uts <- rep_len(units, length(nm)); names(uts) <- nm
  key <- x[[by]]
  FLCore::FLQuants(stats::setNames(lapply(nm, function(.k)
    long2flq(x[key == .k], quant = quant, units = uts[[.k]], years = years)), nm))
}

#' Convert FLQuants back to a long table
#'
#' Inverse of \code{\link{long2flqs}()}: one \code{data.table} with the FLQuant dimension columns, \code{data} and a
#' column \code{name} holding the element names.
#' @param x FLQuants (or FLQuant)
#' @param na.rm drop rows with \code{NA} data
#' @return data.table
#' @export
flqs2long <- function(x, na.rm = TRUE){
  if(is(x, "FLQuant")) x <- FLCore::FLQuants(x = x)
  out <- rbindlist(lapply(names(x), function(.k){
    d <- as.data.table(methods::as(x[[.k]], "data.frame"))
    d[, name := .k]
  }), use.names = TRUE, fill = TRUE)
  if(na.rm) out <- out[!is.na(data)]
  out[]
}

# global variables used in data.table expressions
utils::globalVariables(c(".", "name", "fleet", "year", "unit", "season", "area", "iter", "data", "len", "age",
  "obs", "exp", "se", "se_input", "nsamp", "nsamp_adj", "nsamp_in", "part", "sexes", "sex", "bin", "ageerr",
  "use", "month", "pearson", "type", "i.fleet", "i.name", "i.part", "i.units", "Yr", "Seas", "Fleet", "n", "w", "partname", "qname", "nexp", "np", "ne",
  "used", "effN", "units", "se_total", "nlen", "time", "dev", "q"))

# an r4ss::SS_output() list
.is_ss3rep <- function(x) is.list(x) && !is.data.frame(x) && !is.null(x$FleetNames)

# effective end year: retrospective runs keep endyr and report the peel end year in Retro_year
.endyr <- function(ss3rep){
  ry <- ss3rep$Retro_year
  if(!is.null(ry) && is.finite(ry) && ry < ss3rep$endyr) ry else ss3rep$endyr
}
