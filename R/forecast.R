# forecast.R - short-term forecast scenarios run with SS3 (forecast.ss with fixed F or catch by fleet)

#' Forecast scenarios from SS3 runs
#'
#' Builds FLStockR objects for a set of SS3 forecast runs, e.g. F scenarios (Fadv, Fsq, F0) set up in
#' \code{forecast.ss} and read with \code{r4ss::SSgetoutput()}: each run is converted with \code{\link{ss3stk}()}
#' over the assessment and forecast years, and the agreed benchmarks can replace the SS3 reference points.
#'
#' @param x named list of SS3 reports (one per scenario)
#' @param refpts optional FLPar of benchmarks (e.g. \code{Ftgt}, \code{Btgt}, \code{Bpa}, \code{Blim},
#'   \code{Btrigger}) replacing the SS3 reference points (labels then Btgt / Ftgt)
#' @param output \code{"mle"} (default) or \code{"iters"} (MVLN draws, \code{iters})
#' @param iters number of draws
#' @param years years (default from the start year to the last forecast year with SSB, F and recruitment in the
#'   derived quantities)
#' @param ... further arguments to \code{\link{ss3stk}()} (e.g. \code{Fref})
#' @return FLStocks named by scenario, with the assessment end year as attribute \code{"endyr"}
#' @seealso \code{\link{ss3plot_fwd}()}, \code{\link{ss3_fwdtable}()}
#' @examples
#' fw <- readRDS(system.file("extdata", "hke.37.12-16", "hke.37.12-16_fwd.rds", package = "ss3diags2"))
#' fs <- ss3_fwd(fw)
#' ss3plot_fwd(fs)
#' @export
ss3_fwd <- function(x, refpts = NULL, output = c("mle", "iters"), iters = 1000, years = NULL, ...){
  output <- match.arg(output)
  if(.is_ss3rep(x)) x <- list(x)
  nm <- names(x); if(is.null(nm) || any(nm == "")) nm <- paste0("F", seq_along(x)); names(x) <- nm
  out <- lapply(nm, function(n){
    r <- x[[n]]
    lab <- r$derived_quants$Label
    fy <- function(p) suppressWarnings(as.numeric(sub(p, "", grep(paste0("^", p, "[0-9]{4}$"), lab, value = TRUE))))
    yrs <- years %||% seq(r$startyr, min(max(fy("SSB_")), max(fy("F_")), max(fy("Recr_"))))
    s <- ss3stk(r, output = output, iters = iters, years = yrs, name = n, ...)
    attr(s, "endyr") <- r$endyr
    if(!is.null(refpts)){
      FLCore::refpts(s) <- refpts
      attr(s, "basis") <- utils::modifyList(attr(s, "basis") %||% list(), list(B = "Btgt", F = "Ftgt", C = "MSY"))
      attr(s, "status") <- NULL; attr(s, "mle") <- NULL
    }
    s
  })
  out <- FLCore::FLStocks(stats::setNames(out, nm))
  attr(out, "endyr") <- x[[1]]$endyr
  out
}

# number of intermediate years: forecast years before the scenarios differ in F (0 with one scenario)
.fwd_nint <- function(x, ey){
  if(length(x) < 2) return(0)
  yrs <- as.numeric(dimnames(x[[1]])$year)
  # MLE F where available (MVLN draws differ between runs), else medians
  f <- sapply(as.list(x), function(s){ m <- attr(s, "mle")
    if(is.data.frame(m) && "F" %in% names(m)) m$F[match(yrs, m$year)] else as.numeric(FLCore::iterMedians(FLCore::fbar(s))) })
  d <- which(yrs > ey & apply(f, 1, function(v) diff(range(v, na.rm = TRUE)) > 1e-3*max(abs(v), na.rm = TRUE)))
  if(length(d)) yrs[d[1]] - ey - 1 else 0
}

#' Plot forecast scenarios
#'
#' Short-term forecast plot as in the GFCM / ICES advice: the assessment (black) up to the end year and the forecast
#' scenarios (colours) from the end year on, with or without uncertainty bands, the reference points and vertical
#' lines at the end of the assessment and the start of the advice year.
#'
#' @param x FLStocks from \code{\link{ss3_fwd}()} (or a named list of SS3 forecast reports)
#' @param stk optional FLStockR of the assessment (default: the first scenario up to the end year)
#' @param type quantities, as in \code{\link{ss3plot_trj}()} (default recruitment, SSB, F, catch)
#' @param start first year shown (default 10 years before the end year)
#' @param end last year shown (default all forecast years)
#' @param ci draw uncertainty bands (needs iterations, \code{ss3_fwd(output = "iters")}, or MLE standard errors)
#' @param vlines years after the end year for the vertical dashed lines (default: the end year and the last
#'   intermediate year, i.e. before the first year in which the scenarios differ)
#' @param colours scenario colours (default r4ss colours)
#' @param probs quantiles for the bands
#' @param ncol facet columns
#' @param ... further arguments to \code{\link{ss3plot_trj}()}
#' @return ggplot
#' @export
ss3plot_fwd <- function(x, stk = NULL, type = c("Rec", "SSB", "F", "Catch"), start = NULL, end = NULL, ci = TRUE,
                        vlines = NULL, colours = NULL, probs = c(0.05, 0.2, 0.5, 0.8, 0.95), ncol = 2, ...){
  if(is.list(x) && !is(x, "FLStocks") && .is_ss3rep(x[[1]])) x <- ss3_fwd(x)
  ey <- attr(x, "endyr") %||% attr(x[[1]], "endyr")
  start <- start %||% (ey - 10)
  vlines <- vlines %||% unique(c(0, .fwd_nint(x, ey)))
  a <- .win_stk(stk %||% x[[1]], start = start, end = ey)
  attr(a, "endyr") <- ey
  sc <- lapply(as.list(x), function(s){ w <- .win_stk(s, start = ey, end = end); attr(w, "endyr") <- ey; w })
  st <- FLCore::FLStocks(c(list(Assessment = a), sc))
  cols <- c(Assessment = "black", stats::setNames(colours %||% ss3runcols(length(sc)), names(sc)))
  p <- ss3plot_trj(st, type = type, probs = probs, colours = cols, ncol = ncol, fwd = TRUE, ci = ci, ...)
  # replace the generic forecast marker by the advice lines
  p$layers <- Filter(function(l) !inherits(l$geom, "GeomVline"), p$layers)
  p$labels$caption <- NULL
  if(length(vlines)) p <- p + ggplot2::geom_vline(xintercept = ey + vlines, linetype = 2, linewidth = 0.3)
  p + ggplot2::labs(colour = NULL)
}

#' Fishing opportunities table of forecast scenarios
#'
#' Summary of forecast scenarios as in the GFCM / ICES advice tables: catch and F in the advice year, SSB at the
#' start of the following year, its change relative to the advice year, and the probability that it is below Blim
#' (from iterations).
#'
#' @param x FLStocks from \code{\link{ss3_fwd}()}
#' @param year advice year (default: the first year in which the scenarios differ in F, after the intermediate
#'   years)
#' @param unc optional FLStocks with iterations for the probabilities (default \code{x} when it has iterations)
#' @param blim Blim (default the \code{Blim} reference point of the stocks)
#' @return data.table with Basis, catch, F, SSB(year + 1), SSB change (\%) and P(SSB < Blim)
#' @export
ss3_fwdtable <- function(x, year = NULL, unc = NULL, blim = NULL){
  ey <- attr(x, "endyr") %||% attr(x[[1]], "endyr")
  y <- as.character(year %||% (ey + .fwd_nint(x, ey) + 1)); y1 <- as.character(as.numeric(y) + 1)
  med <- function(q) as.numeric(FLCore::iterMedians(q))
  unc <- unc %||% (if(FLCore::dims(x[[1]])$iter > 1) x else NULL)
  rbindlist(lapply(names(x), function(n){
    s <- x[[n]]
    if(!y1 %in% dimnames(s)$year) stop("year ", y1, " not in the forecast of ", n, call. = FALSE)
    b1 <- med(FLCore::ssb(s)[, y1]); b0 <- med(FLCore::ssb(s)[, y])
    bl <- blim %||% .rp(s, "Blim")
    pb <- if(!is.null(unc) && is.finite(bl)) mean(as.numeric(FLCore::ssb(unc[[n]])[, y1]) < bl) else NA_real_
    d <- data.table(Basis = n, C = round(med(FLCore::catch(s)[, y]), 1), F = round(med(FLCore::fbar(s)[, y]), 3),
                    SSB = round(b1), SSB_change = round(100*(b1/b0 - 1), 1), P_Blim = round(pb, 3))
    setnames(d, c("C", "F", "SSB", "P_Blim"), c(paste0("C_", y), paste0("F_", y), paste0("SSB_", y1),
                                                 paste0("P(SSB_", y1, "<Blim)")))
    d
  }))
}
