# Ensemble catch advice from a grid of catch (TAC) forecasts

#' Summary of a catch scenario in the advice year
#'
#' Ensemble (or single model) medians and probabilities of one forecast scenario: catch and F/Ftgt in the advice
#' year, SSB at the start of the advice year and of the next year, and the probabilities of SSB below Blim and
#' Btrigger in the next year. Ratios and probabilities use each draw's own reference points (each model's own for an
#' ensemble).
#'
#' @param s FLStockR with iterations (one scenario of \code{\link{ss3_fwd}()} or \code{\link{ss3_fwdensemble}()})
#' @param y advice year
#' @param btrigger Btrigger as a fraction of Btgt
#' @return one-row data.table
#' @noRd
.scen_summary <- function(s, y, btrigger = 1){
  y <- as.character(y); y1 <- as.character(as.numeric(y) + 1)
  b0 <- as.numeric(FLCore::ssb(s)[, y]); b1 <- as.numeric(FLCore::ssb(s)[, y1])
  rp <- FLCore::refpts(s); n <- length(b1)
  rpv <- function(p) if(p %in% dimnames(rp)$params) rep_len(as.numeric(rp[p, ]), n) else rep(NA_real_, n)
  bt <- rpv("Btgt"); bl <- rpv("Blim"); btr <- btrigger*bt
  st <- attr(s, "status")
  fr <- if(!is.null(st$harvest)) as.numeric(st$harvest[, y]) else as.numeric(FLCore::fbar(s)[, y])/rpv("Ftgt")
  # ratios to the targets from the ratio draws (within each model, as SSdeltaMVLN), else from SSB and the benchmarks
  r0 <- if(!is.null(st$stock)) as.numeric(st$stock[, y]) else b0/bt
  r1 <- if(!is.null(st$stock)) as.numeric(st$stock[, y1]) else b1/bt
  data.table(catch = stats::median(as.numeric(FLCore::catch(s)[, y])), F = stats::median(as.numeric(FLCore::fbar(s)[, y])),
             Fratio = stats::median(fr), SSB0 = stats::median(b0), SSB1 = stats::median(b1),
             B0ratio = stats::median(r0)/btrigger, B1ratio = stats::median(r1)/btrigger,
             pBlim = mean(r1 < bl/bt), pBtrigger = mean(r1 < btrigger))
}

#' Value at which a monotone quantity of the scenarios reaches a target
#'
#' @param x scenario catches (increasing)
#' @param y quantity by scenario
#' @param target target value
#' @return interpolated catch at the first crossing, or NA when the target is not bracketed
#' @noRd
.cross <- function(x, y, target){
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  j <- which(diff(sign(y - target)) != 0 | (y[-length(y)] - target) == 0)
  if(!length(j)) return(NA_real_)
  j <- j[1]
  if(y[j] == target) return(x[j])
  x[j] + (target - y[j])/(y[j + 1] - y[j])*(x[j + 1] - x[j])
}

#' Catch options from a grid of catch scenarios
#'
#' Catch options for the advice year as in the ICES catch scenario tables, from the forecasts of a fine grid of catch
#' (TAC) scenarios of one model or a model ensemble (\code{\link{ss3_fwdensemble}()}). Each option is the catch at
#' which an ensemble median (or probability) reaches its target, interpolated linearly between neighbouring
#' scenarios; the other quantities are interpolated at that catch. Ratios to the reference points and probabilities
#' are computed within each model (each draw's own reference points), so F/Fmsy = 1 is the "average Fmsy" of the
#' ensemble.
#'
#' The options (\code{bases}):
#' \describe{
#'   \item{\code{"Fmsy"}, \code{"Flower"}, \code{"Fupper"}}{F/Ftgt = 1, \code{flower}, \code{fupper}, reduced by
#'     SSB(advice year)/Btrigger when below 1 (\code{rule = TRUE}: ICES advice rule and EU MAP)}
#'   \item{\code{"F0"}}{the lowest catch scenario}
#'   \item{\code{"Fsq"}}{F of the advice year equal to F of the last intermediate year}
#'   \item{\code{"Fpa"}}{F/Ftgt = \code{fpa}}
#'   \item{\code{"SSB=SSB"}}{SSB at the start of the next year equal to SSB at the start of the advice year}
#'   \item{\code{"SSB=Btrigger"}, \code{"SSB=Bpa"}, \code{"SSB=Blim"}}{median SSB of the next year at the reference
#'     point}
#'   \item{\code{"P(Blim)"}}{P(SSB next year < Blim) = \code{risk}}
#' }
#'
#' @param x FLStocks of catch scenarios with iterations: \code{\link{ss3_fwdensemble}()} or
#'   \code{ss3_fwd(output = "iters")}
#' @param year advice year (default: the first year in which the scenarios differ)
#' @param bases catch options (see Details)
#' @param btrigger MSY Btrigger as a fraction of Btgt (default 1: Btrigger = Btgt, e.g. B30\%)
#' @param bpa Bpa as a fraction of Btgt (default \code{btrigger})
#' @param flower,fupper,fpa F/Fmsy of the F range and of Fpa
#' @param rule apply the advice rule (F reduced by SSB/Btrigger below Btrigger) to the Fmsy, Flower and Fupper options
#' @param risk target risk for \code{"P(Blim)"}
#' @param advice previous advice (catch) for the \% advice change: one value, or named by basis (e.g. the EU MAP
#'   advice for Fmsy, Flower and Fupper)
#' @param tac catch (TAC) of the intermediate year for the \% TAC change (default: the intermediate year's catch of
#'   the scenarios)
#' @return data.table: basis, catch, F/Fmsy, SSB (next year), SSB/Btrigger, \% SSB change, \% TAC change, \% advice
#'   change, P(SSB < Blim) and P(SSB < Btrigger) (\%), with attributes \code{year} and \code{scenarios} (the summary
#'   of every scenario)
#' @seealso \code{\link{ss3plot_catchopts}()}, \code{\link{ss3_fwdprobs}()}, \code{\link{ss3_fwdrisk}()}
#' @examples
#' \donttest{
#' fw <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18_fwd.rds", package = "ss3diags2"))
#' fs <- ss3_fwd(fw, output = "iters", iters = 500, Fref = "Btgt", blim = 0.25, bpa = 0.5)
#' ss3_catchopts(fs, bases = c("F0", "Fmsy", "Fsq"), btrigger = 0.5)
#' }
#' @export
ss3_catchopts <- function(x, year = NULL, bases = c("Fmsy", "Flower", "Fupper", "F0", "Fsq", "SSB=SSB",
                                                    "SSB=Btrigger", "SSB=Bpa", "SSB=Blim", "P(Blim)"),
                          btrigger = 1, bpa = btrigger, flower = 0.75, fupper = 1.21, fpa = NULL, rule = TRUE,
                          risk = 0.05, advice = NULL, tac = NULL){
  ey <- attr(x, "endyr") %||% attr(x[[1]], "endyr")
  y <- year %||% (ey + .fwd_nint(x, ey) + 1); y1 <- y + 1
  sm <- rbindlist(lapply(as.list(x), .scen_summary, y = y, btrigger = btrigger))
  sm[, scenario := names(x)]
  setorder(sm, catch)
  # intermediate year: F and catch (same in all scenarios)
  yi <- y - 1
  fi <- stats::median(as.numeric(FLCore::fbar(x[[1]])[, as.character(yi)]))
  ftg <- .rp(x[[1]], "Ftgt")
  ci <- tac %||% stats::median(as.numeric(FLCore::catch(x[[1]])[, as.character(yi)]))
  b0r <- sm$B0ratio[1]                                          # SSB(advice year)/Btrigger, the same in all scenarios
  red <- if(rule) min(1, b0r) else 1
  bl <- .rp(x[[1]], "Blim")/.rp(x[[1]], "Btgt")/btrigger       # Blim as fraction of Btrigger
  tg <- list(Fmsy = c("Fratio", red), Flower = c("Fratio", flower*red), Fupper = c("Fratio", fupper*red),
             Fpa = c("Fratio", fpa %||% NA), Fsq = c("F", fi), `SSB=SSB` = c("SSBchange", 0),
             `SSB=Btrigger` = c("B1ratio", 1), `SSB=Bpa` = c("B1ratio", bpa/btrigger), `SSB=Blim` = c("B1ratio", bl),
             `P(Blim)` = c("pBlim", risk))
  sm[, SSBchange := SSB1/SSB0 - 1]
  cols <- c("F", "Fratio", "SSB1", "B1ratio", "SSBchange", "pBlim", "pBtrigger")
  at <- function(cc) { if(!is.finite(cc)) return(as.list(stats::setNames(rep(NA_real_, length(cols)), cols)))
                       lapply(stats::setNames(cols, cols), function(k) stats::approx(sm$catch, sm[[k]], cc, ties = mean)$y) }
  out <- rbindlist(lapply(bases, function(b){
    if(b == "F0"){ cc <- sm$catch[1] } else {
      if(!b %in% names(tg)) stop("unknown basis '", b, "'", call. = FALSE)
      k <- tg[[b]][1]; v <- as.numeric(tg[[b]][2])
      cc <- if(is.finite(v)) .cross(sm$catch, sm[[k]], v) else NA_real_
    }
    data.table(basis = b, catch = cc, as.data.table(at(cc)))
  }))
  adv <- if(is.null(advice)) NA_real_ else if(is.null(names(advice))) advice else advice[out$basis]
  res <- out[, .(basis, catch = round(catch), `F/Fmsy` = round(Fratio, 2), SSB = round(SSB1),
                 `SSB/Btrigger` = round(B1ratio, 2), `SSB change (%)` = round(100*SSBchange, 1),
                 `TAC change (%)` = round(100*(catch/ci - 1), 1), `advice change (%)` = round(100*(catch/adv - 1), 1),
                 `P(SSB<Blim) (%)` = round(100*pBlim, 2), `P(SSB<Btrigger) (%)` = round(100*pBtrigger, 1))]
  setnames(res, c("catch", "F/Fmsy", "SSB", "SSB/Btrigger", "P(SSB<Blim) (%)", "P(SSB<Btrigger) (%)"),
           c(paste0("catch ", y), paste0("F", y, "/Fmsy"), paste0("SSB", y1), paste0("SSB", y1, "/Btrigger"),
             paste0("P(SSB", y1, "<Blim) (%)"), paste0("P(SSB", y1, "<Btrigger) (%)")))
  if(all(is.na(adv))) res[, `advice change (%)` := NULL]
  attr(res, "year") <- y
  attr(res, "scenarios") <- sm[]
  attr(res, "interim") <- data.table(year = yi, catch = ci, F = fi, `F/Fmsy` = fi/ftg, `SSB/Btrigger` = b0r)
  res[]
}

#' Plot catch options
#'
#' F/Fmsy, SSB/Btrigger in the year after the advice year and P(SSB < Blim) against the catch of the scenarios (ensemble
#' medians), with the catch options of \code{\link{ss3_catchopts}()} marked: the harvest strategy view of a catch
#' scenario grid.
#'
#' @param x output of \code{\link{ss3_catchopts}()}
#' @param bases options to mark (default all with a catch)
#' @param risk target risk line
#' @param text show the title, subtitle and caption (\code{FALSE}: none; labels inside the plot and legends stay)
#' @return ggplot
#' @examples
#' \donttest{
#' fw <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18_fwd.rds", package = "ss3diags2"))
#' fs <- ss3_fwd(fw, output = "iters", iters = 500, Fref = "Btgt", blim = 0.25, bpa = 0.5)
#' ss3plot_catchopts(ss3_catchopts(fs, bases = c("F0", "Fmsy", "Fsq"), btrigger = 0.5))
#' }
#' @export
ss3plot_catchopts <- function(x, bases = NULL, risk = 0.05, text = TRUE){
  y <- attr(x, "year"); sm <- copy(attr(x, "scenarios"))
  d <- rbind(sm[, .(catch, panel = paste0("F", y, "/Fmsy"), value = Fratio, ref = 1)],
             sm[, .(catch, panel = paste0("SSB", y + 1, "/Btrigger"), value = B1ratio, ref = 1)],
             sm[, .(catch, panel = paste0("P(SSB", y + 1, " < Blim)"), value = pBlim, ref = risk)])
  d[, panel := factor(panel, levels = unique(panel))]
  o <- x[is.finite(get(paste0("catch ", y)))]
  if(!is.null(bases)) o <- o[basis %in% bases]
  o <- data.table(basis = o$basis, catch = o[[paste0("catch ", y)]])
  oc <- stats::setNames(ss3cols(max(nrow(o), 2))[seq_len(nrow(o))], o$basis)
  p <- ggplot2::ggplot(d, ggplot2::aes(catch, value)) +
    ggplot2::geom_hline(ggplot2::aes(yintercept = ref), linetype = 2, colour = "grey40") +
    ggplot2::geom_line(linewidth = 0.9) + ggplot2::geom_point(size = 0.8) +
    ggplot2::geom_vline(data = o, ggplot2::aes(xintercept = catch, colour = basis), linewidth = 0.6) +
    ggplot2::scale_colour_manual(values = oc, name = NULL) +
    ggplot2::facet_wrap(~ panel, ncol = 1, scales = "free_y") +
    ggplot2::labs(x = paste("Catch", y), y = NULL, title = "Catch options",
                  subtitle = "ensemble medians and probability by catch scenario; vertical lines: catch options") +
    theme_ss3()
  .drop_text(p, text)
}

#' Probability tables of catch scenarios
#'
#' Probabilities by catch scenario and year, as the Kobe strategy matrices of ensemble advice (\code{kobe:::smry}):
#' P(F < Ftgt), P(SSB > Btgt), P(green: SSB > Btgt and F < Ftgt), P(SSB > Btrigger) and P(SSB > Blim), computed with
#' each draw's own reference points.
#'
#' @param x FLStocks of catch scenarios with iterations (\code{\link{ss3_fwdensemble}()} or
#'   \code{ss3_fwd(output = "iters")})
#' @param years years (default: the intermediate years, the advice year and the next year)
#' @param btrigger Btrigger as a fraction of Btgt
#' @param long return a long table instead of one wide table (scenario by year) per probability
#' @return named list of data.tables (scenario x year, in \%), or a long data.table
#' @examples
#' \donttest{
#' fw <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18_fwd.rds", package = "ss3diags2"))
#' fs <- ss3_fwd(fw, output = "iters", iters = 500, Fref = "Btgt", blim = 0.25, bpa = 0.5)
#' ss3_fwdprobs(fs, btrigger = 0.5)
#' }
#' @export
ss3_fwdprobs <- function(x, years = NULL, btrigger = 1, long = FALSE){
  ey <- attr(x, "endyr") %||% attr(x[[1]], "endyr")
  years <- years %||% seq(ey + 1, ey + .fwd_nint(x, ey) + 2)
  d <- rbindlist(lapply(names(x), function(n){
    s <- x[[n]]; rp <- FLCore::refpts(s); st <- attr(s, "status")
    rbindlist(lapply(as.character(years[as.character(years) %in% dimnames(s)$year]), function(y){
      b <- as.numeric(FLCore::ssb(s)[, y]); nn <- length(b)
      rpv <- function(p) if(p %in% dimnames(rp)$params) rep_len(as.numeric(rp[p, ]), nn) else rep(NA_real_, nn)
      fr <- if(!is.null(st$harvest)) as.numeric(st$harvest[, y]) else as.numeric(FLCore::fbar(s)[, y])/rpv("Ftgt")
      br <- if(!is.null(st$stock)) as.numeric(st$stock[, y]) else b/rpv("Btgt")
      data.table(scenario = n, year = as.numeric(y), `F<Ftgt` = mean(fr < 1), `SSB>Btgt` = mean(br > 1),
                 green = mean(fr < 1 & br > 1), `SSB>Btrigger` = mean(br > btrigger),
                 `SSB>Blim` = mean(br > rpv("Blim")/rpv("Btgt")))
    }))
  }))
  if(long) return(d[])
  pr <- setdiff(names(d), c("scenario", "year"))
  stats::setNames(lapply(pr, function(k){
    w <- dcast(d, scenario ~ year, value.var = k)
    w <- w[match(names(x), scenario)]
    for(j in setdiff(names(w), "scenario")) set(w, j = j, value = round(100*w[[j]], 1))
    w }), pr)
}

utils::globalVariables(c("Fratio", "B1ratio", "SSB0", "SSB1", "SSBchange", "pBlim", "pBtrigger", "panel", "ref",
                         "basis"))
