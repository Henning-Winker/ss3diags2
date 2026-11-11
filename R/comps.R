# comps.R - length and age compositions by fleet

.part_lab <- c("0" = "all", "1" = "discard", "2" = "retained")

# long table of one SS3 composition database (lendbase, agedbase)
.ss3_comps <- function(ss3rep, db, quant, fleets, what, scale, nsamp, sex, part, long){
  db <- as.data.table(db)
  if(!nrow(db)) stop("no ", quant, " compositions in this SS3 report", call. = FALSE)
  # r4ss splits marginal age compositions (agedbase) from conditional age-at-length (condbase)
  fleets <- .pick_fleets(ss3rep, fleets, .fleetnames(ss3rep)[sort(unique(db$Fleet))])
  ns <- ss3rep$nseasons %||% 1; na <- ss3rep$nareas %||% 1; nsx <- ss3rep$nsexes %||% 1
  nin <- .col(db, c("Nsamp_in", "N"))
  if(nsamp == "input" && is.null(nin)){
    message("no input sample size in this r4ss output; using Nsamp_adj"); nsamp <- "adj"
  }
  d <- data.table(fleet = db$Fleet, name = .fleetnames(ss3rep)[db$Fleet], bin = db$Bin, year = db$Yr,
                  unit = .unit_lab(.col(db, "Sex", 1), .col(db, "Sexes", 0), nsx),
                  season = .season_lab(db$Seas, ns), area = .area_lab(.col(db, "Area", 1), na), iter = "1",
                  part = .col(db, "Part", 0), ageerr = .col(db, "Ageerr", 0), month = .col(db, "Month", NA_real_),
                  obs = db$Obs, exp = db$Exp, pearson = .col(db, "Pearson", NA_real_),
                  nsamp_adj = db$Nsamp_adj, nsamp_in = nin %||% NA_real_, effN = .col(db, "effN", NA_real_),
                  used = .col(db, "Used", "yes"))
  d <- d[name %in% fleets]
  if(!is.null(part)){ keep <- part; d <- d[part %in% keep] }
  d[, nsamp := if(nsamp == "adj") nsamp_adj else nsamp_in]
  # one FLQuant per fleet and partition (ageing-error definitions, often year-specific, stay in column ageerr)
  d[, partname := name]
  d[, np := uniqueN(part), by = name][np > 1, partname := paste0(name, ".", .part_lab[as.character(part)])]
  d[, np := NULL]
  d[, n := obs*nsamp]
  d[, nexp := exp*nsamp]
  d[, data := switch(what, obs = if(scale == "n") n else obs, exp = if(scale == "n") nexp else exp,
                           pearson = pearson)]
  setnames(d, "bin", quant)
  by <- c("partname", quant, .flqdims)
  # several samples in one cell (sexes combined, months within a season, replicates): numbers are summed and
  # proportions recomputed from the summed numbers within each composition (partname x year x unit x season x area)
  if(sex == "combine") d[, unit := "unique"]
  if(anyDuplicated(d, by = by)){
    d <- d[, .(fleet = fleet[1], name = name[1], part = part[1], ageerr = ageerr[1], month = mean(month),
               pearson = mean(pearson), n = sum(n), nexp = sum(nexp), nsamp_adj = mean(nsamp_adj),
               nsamp_in = mean(nsamp_in), effN = mean(effN), used = used[1]), by = by]
    d[, `:=`(obs = n/sum(n, na.rm = TRUE), exp = nexp/sum(nexp, na.rm = TRUE)), by = c("partname", .flqdims)]
    d[, data := switch(what, obs = if(scale == "n") n else obs, exp = if(scale == "n") nexp else exp,
                             pearson = pearson)]
  }
  d[, units := if(what == "pearson") "NA" else if(scale == "n") "numbers" else "prop"]
  setcolorder(d, c("fleet", "name", "partname", quant, .flqdims, "data"))
  if(long) return(d[])
  nms <- unique(d$partname)
  long2flqs(d, quant = quant, units = d$units[1], by = "partname", order = nms)
}

#' Length and age compositions by fleet
#'
#' Extracts SS3 length (\code{ss3rep$lendbase}) or marginal age (\code{ss3rep$agedbase}) compositions as
#' \code{FLQuants}, one \code{FLQuant} per fleet, with length bins (lower bounds, as in SS3) or ages as the first
#' dimension, sex as \code{unit} ("F", "M"; "unique" for combined-sex samples or single-sex models), season and
#' area as in the model. Years without samples are \code{NA}. Length compositions by fleet are the \code{lfd} input
#' of FLicc (see \code{\link{ss2flicc}()} for the exact shape).
#'
#' SS3 reports compositions as proportions with a sample size; \code{scale = "n"} multiplies by the sample size
#' (\code{nsamp = "adj"}: after variance adjustment, \code{"input"}: as input) to give numbers. Fleets with more than
#' one partition (discards, retained) get one element each, named \code{fleet.discard}, \code{fleet.retained}.
#' Ageing-error definitions (often year-specific) are kept in the long table (column \code{ageerr}).
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param fleets fleet numbers or names (default: all fleets with compositions)
#' @param what \code{"obs"} (default), \code{"exp"} (expected) or \code{"pearson"} (Pearson residuals)
#' @param scale \code{"n"} (default, numbers = proportion x sample size) or \code{"prop"} (proportions)
#' @param nsamp \code{"adj"} (default) or \code{"input"} sample size for \code{scale = "n"}
#' @param sex \code{"unit"} (default: sexes in the unit dimension) or \code{"combine"} (summed over sexes)
#' @param part optional SS3 partitions to keep (0 = all, 1 = discard, 2 = retained)
#' @param long return the long data.table instead of FLQuants
#' @return FLQuants (named by fleet), or a data.table with the FLQuant dimension columns (\code{len} or \code{age},
#'   \code{year}, \code{unit}, \code{season}, \code{area}, \code{iter}), \code{data} and the SS3 columns
#'   \code{fleet}, \code{name}, \code{partname}, \code{part}, \code{ageerr}, \code{month}, \code{obs}, \code{exp},
#'   \code{pearson}, \code{nsamp}, \code{nsamp_adj}, \code{nsamp_in}, \code{effN}, \code{n}, \code{units}
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' lfd <- ss3_lencomp(ss3sma)                    # numbers at length, sexes as unit
#' lfd <- ss3_lencomp(ss3sma, sex = "combine")   # sexes combined
#' data(pac.hke, package = "ss3diags")
#' afd <- ss3_agecomp(ss3phk, scale = "prop")
#' }
#' @export
ss3_lencomp <- function(ss3rep, fleets = NULL, what = c("obs", "exp", "pearson"), scale = c("n", "prop"),
                        nsamp = c("adj", "input"), sex = c("unit", "combine"), part = NULL, long = FALSE){
  .ss3_comps(ss3rep, ss3rep$lendbase, "len", fleets, match.arg(what), match.arg(scale), match.arg(nsamp),
             match.arg(sex), part, long)
}

#' @rdname ss3_lencomp
#' @export
ss3_agecomp <- function(ss3rep, fleets = NULL, what = c("obs", "exp", "pearson"), scale = c("n", "prop"),
                        nsamp = c("adj", "input"), sex = c("unit", "combine"), part = NULL, long = FALSE){
  .ss3_comps(ss3rep, ss3rep$agedbase, "age", fleets, match.arg(what), match.arg(scale), match.arg(nsamp),
             match.arg(sex), part, long)
}
