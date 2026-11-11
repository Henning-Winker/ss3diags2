# indices.R - abundance indices as FLIndices

#' Abundance indices
#'
#' Extracts the abundance indices fitted in SS3 (\code{ss3rep$cpue}) as \code{FLIndices} of \code{FLIndexBiomass},
#' one per fleet, the index input of JABBA3. Indices in numbers and in weight are both stored as
#' \code{FLIndexBiomass}; the SS3 index units are recorded in \code{units(index(x))} ("numbers", "biomass", ...).
#'
#' Slots filled: \code{index} (observed, or expected with \code{what = "exp"}), \code{index.var} (log-scale SE),
#' \code{index.q} (SS3 \code{Calc_Q}), \code{range} \code{startf}/\code{endf} (timing within the year from the SS3
#' month), \code{distribution} ("lognormal"), \code{name} and \code{desc}. With \code{long = TRUE} the full table
#' is returned (observed, expected, SEs, deviations, catchability, use flag), the basis for residual diagnostics.
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param fleets fleet numbers or names (default: all fleets with index data)
#' @param what \code{"obs"} (default) or \code{"exp"} (fitted index) in the \code{index} slot
#' @param se \code{"input"} (default): the input SE (\code{SE_input}), as JABBA3 adds its own estimated SD
#'   (\code{sdi}) to it, or \code{"total"}: the SE used in the SS3 likelihood (input SE plus SS3 extra SD).
#'   Older r4ss outputs without \code{SE_input} only have the total SE.
#' @param excluded \code{"drop"} (default): observations excluded in SS3 (negative fleet, \code{Use} < 0) are set to
#'   NA; \code{"keep"}: kept
#' @param long return the long data.table instead of FLIndices
#' @return FLIndices (named by fleet), or a data.table with columns \code{fleet}, \code{name}, \code{age} ("all"),
#'   \code{year}, \code{unit}, \code{season}, \code{area}, \code{iter}, \code{data} and \code{month}, \code{time},
#'   \code{obs}, \code{exp}, \code{se}, \code{se_input}, \code{se_total}, \code{q}, \code{dev}, \code{use},
#'   \code{units}
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' idx <- ss3_indices(ss3sma)
#' names(idx)
#' index.var(idx[[1]])
#' }
#' @export
ss3_indices <- function(ss3rep, fleets = NULL, what = c("obs", "exp"), se = c("input", "total"),
                        excluded = c("drop", "keep"), long = FALSE){
  what <- match.arg(what); se <- match.arg(se); excluded <- match.arg(excluded)
  cp <- as.data.table(ss3rep$cpue)
  if(!nrow(cp)) stop("no index data in ss3rep$cpue", call. = FALSE)
  fl <- ss3_fleets(ss3rep)
  fleets <- .pick_fleets(ss3rep, fleets, .fleetnames(ss3rep)[sort(unique(cp$Fleet))])
  ns <- ss3rep$nseasons %||% 1; na <- ss3rep$nareas %||% 1
  se_in <- .col(cp, "SE_input")
  if(se == "input" && is.null(se_in))
    message("ss3_indices: no SE_input in this r4ss output; index.var holds the total SE (input + extra SD)")
  mon <- .col(cp, "Month", NA_real_)
  d <- data.table(fleet = cp$Fleet, name = .fleetnames(ss3rep)[cp$Fleet], age = "all", year = cp$Yr, unit = "unique",
                  season = .season_lab(cp$Seas, ns), area = .area_lab(.col(cp, "Area", 1), na), iter = "1",
                  month = mon, time = .col(cp, "Time", cp$Yr), obs = cp$Obs, exp = cp$Exp,
                  se_input = se_in %||% NA_real_, se_total = cp$SE, q = .col(cp, "Calc_Q", NA_real_),
                  dev = .col(cp, "Dev", NA_real_), use = .col(cp, "Use", 1))
  d <- d[name %in% fleets]
  d[, se := if(se == "input" && !is.null(se_in)) se_input else se_total]
  if(excluded == "drop") d <- d[is.na(use) | use > 0]
  d[, data := if(what == "obs") obs else exp]
  d[, units := fl$survey_units[match(fleet, fl$fleet)]]
  setcolorder(d, c("fleet", "name", "age", "year", "unit", "season", "area", "iter", "data"))
  if(anyDuplicated(d, by = c("name", "year", "season", "area")))
    stop("ss3_indices: several observations of one index in the same year/season/area (sub-seasons?); ",
         "use long = TRUE", call. = FALSE)
  if(long) return(d[])
  FLCore::FLIndices(stats::setNames(lapply(fleets, function(f){
    x <- d[name == f]
    idx <- long2flq(x, quant = "age", units = x$units[1])
    v <- long2flq(x[, .(age, year, unit, season, area, iter, data = se)], quant = "age", units = "se",
                  years = dimnames(idx)$year)
    q <- long2flq(x[, .(age, year, unit, season, area, iter, data = q)], quant = "age", units = "NA",
                  years = dimnames(idx)$year)
    fi <- FLCore::FLIndexBiomass(index = idx, name = f,
                                 desc = paste0("SS3 fleet ", x$fleet[1], " (", x$units[1], "), ", what))
    FLCore::index.var(fi) <- v
    FLCore::index.q(fi) <- q
    fi@distribution <- "lognormal"
    tf <- if(all(is.na(x$month))) mean(x$time - x$year) else mean((x$month - 1)/12, na.rm = TRUE)
    fi@range[c("startf", "endf")] <- tf
    fi
  }), fleets))
}
