# forecast.R - short-term forecast scenarios run with SS3 (forecast.ss with fixed F or catch by fleet)

#' Convert \code{r4ss::SSsummarize()} output to minimal SS3 reports
#'
#' R4ss::SSsummarize() output -> list of minimal SS3 reports (derived quantities, years, F basis, targets)
#'
#' @param x \code{SSsummarize()} list
#' @return named list of minimal reports (derived quantities, years, F basis, targets)
#' @noRd
.summ2reps <- function(x, peels = NULL){
  q <- x$quants; qs <- x$quantsSD
  nm <- x$modelnames %||% paste0("model", seq_len(x$n))
  ix <- x$indices
  stats::setNames(lapply(seq_len(x$n), function(i){
    r <- list(derived_quants = data.frame(Label = q$Label, Value = q[[i]], StdDev = qs[[i]]),
              startyr = x$startyrs[i], endyr = x$endyrs[i], F_std_basis = x$FvalueLabels[i] %||% "_abs_F",
              btarg = x$btargs[i], sprtarg = x$sprtargs[i], SpawnOutputUnits = x$SpawnOutputUnits[i],
              nsexes = x$nsexes[i], summary = TRUE)
    fn <- x$FleetNames
    r$FleetNames <- if(is.list(fn)) fn[[i]] else fn
    if(is.data.frame(ix) && "imodel" %in% names(ix)) r$cpue <- ix[ix$imodel == i, , drop = FALSE]
    if(!is.null(peels)) r$Retro_year <- x$endyrs[i] - peels[i]           # retrospective peel end year
    r
  }), nm)
}

#' Retrospective runs as a list of reports
#'
#' A list of SS3 reports, or an \code{r4ss::SSsummarize()} list of the reference run and its peels converted to
#' minimal reports with the peel end years (SS3 retrospective runs keep endyr; the peels are 0, 1, 2, ... years)
#'
#' @param retro list of SS3 reports or \code{SSsummarize()} list (reference run first)
#' @param peels years peeled from each run (default \code{0:(n - 1)})
#' @return list of reports
#' @noRd
.retro_reps <- function(retro, peels = NULL){
  if(.is_summary(retro)) return(.summ2reps(retro, peels %||% (seq_len(retro$n) - 1)))
  if(!is.null(peels)) for(i in seq_along(retro)) retro[[i]]$Retro_year <- retro[[i]]$endyr - peels[i]
  retro
}

#' Is x an \code{r4ss::SSsummarize()} list?
#'
#' @param x object
#' @return logical
#' @noRd
.is_summary <- function(x) is.list(x) && !is.null(x$quants) && !is.null(x$modelnames) && !is.null(x$n)

#' Forecast scenarios from SS3 runs
#'
#' Builds FLStockR objects for a set of SS3 forecast runs, e.g. F scenarios (Fadv, Fsq, F0) set up in
#' \code{forecast.ss} and read with \code{r4ss::SSgetoutput()}: each run is converted with \code{\link{ss3stk}()}
#' over the assessment and forecast years, and the agreed benchmarks can replace the SS3 reference points.
#'
#' @param x named list of SS3 reports (one per scenario), or their summary from \code{r4ss::SSsummarize()} (much
#'   smaller; MLE with standard errors only, and catch for the forecast years only, so pass the assessment to
#'   \code{ss3plot_fwd(stk = )})
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
#' fw <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18_fwd.rds", package = "ss3diags2"))
#' fs <- ss3_fwd(fw, Fref = "Btgt", bpa = 0.5, blim = 0.25)
#' names(fs)
#' hke <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18.rds", package = "ss3diags2"))
#' ss3plot_fwd(fs, stk = ss3stk(hke, Fref = "Btgt", bpa = 0.5, blim = 0.25), ci = FALSE)
#' ss3_fwdtable(fs)
#' @export
ss3_fwd <- function(x, refpts = NULL, output = c("mle", "iters"), iters = 1000, years = NULL, ...){
  # with SSsummarize() input, output = "iters" gives lognormal draws from the MLE and SEs
  output <- match.arg(output)
  draws <- FALSE
  if(.is_summary(x)){
    # SSsummarize() output has no covariance: draws from the MLE and the standard errors (lognormal, see .lnorm_iters)
    draws <- output == "iters"; output <- "mle"; x <- .summ2reps(x)
  }
  if(.is_ss3rep(x) && is.null(x$summary)) x <- list(x)
  nm <- names(x); if(is.null(nm) || any(nm == "")) nm <- paste0("F", seq_along(x)); names(x) <- nm
  out <- lapply(nm, function(n){
    r <- x[[n]]
    lab <- r$derived_quants$Label
    fy <- function(p) suppressWarnings(as.numeric(sub(p, "", grep(paste0("^", p, "[0-9]{4}$"), lab, value = TRUE))))
    yrs <- years %||% seq(r$startyr, min(max(fy("SSB_")), max(fy("F_")), max(fy("Recr_"))))
    s <- ss3stk(r, output = output, iters = iters, years = yrs, name = n, ...)
    attr(s, "endyr") <- r$endyr
    if(draws) s <- .lnorm_iters(s, iters, seed = 1)    # common random numbers: scenarios differ only by their MLE
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

#' Number of intermediate years of forecast scenarios
#'
#' Number of intermediate years: forecast years before the scenarios differ in F (0 with one scenario)
#'
#' @param x FLStocks of scenarios (\code{ss3_fwd()})
#' @param ey assessment end year
#' @return forecast years before the scenarios differ in F (0 with one scenario)
#' @noRd
.fwd_nint <- function(x, ey){
  if(!is.null(attr(x, "nint"))) return(attr(x, "nint"))
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
#' @param vlines years after the end year for the vertical dashed lines (default: the end year and the
#'   implementation year, with the intermediate years between them)
#' @param refnames reference points drawn (default the advice benchmarks Btgt, Bpa, Blim, Btrigger, Ftgt; NULL for
#'   all, including B0, R0 and MSY)
#' @param colours scenario colours (default r4ss colours)
#' @param probs quantiles for the bands
#' @param ncol facet columns
#' @param ... further arguments to \code{\link{ss3plot_trj}()}
#' @param text show the title, subtitle and caption (\code{FALSE}: none; labels inside the plot and legends stay)
#' @return ggplot
#' @examples
#' fw <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18_fwd.rds", package = "ss3diags2"))
#' hke <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18.rds", package = "ss3diags2"))
#' fs <- ss3_fwd(fw, Fref = "Btgt", bpa = 0.5, blim = 0.25)
#' stk <- ss3stk(hke, output = "iters", iters = 200, Fref = "Btgt", bpa = 0.5, blim = 0.25)
#' ss3plot_fwd(fs, stk = stk)                 # with uncertainty
#' ss3plot_fwd(fs, ci = FALSE, start = 2010)  # scenarios only
#' @export
ss3plot_fwd <- function(x, stk = NULL, type = c("Rec", "SSB", "F", "Catch"), start = NULL, end = NULL, ci = TRUE,
                        vlines = NULL, refnames = c("Btgt", "Bpa", "Blim", "Btrigger", "Ftgt"), colours = NULL, probs = c(0.05, 0.2, 0.5, 0.8, 0.95), ncol = 2, ..., text = TRUE){
  .drop_text((function(){
  if(.is_summary(x) || (is.list(x) && !is(x, "FLStocks") && .is_ss3rep(x[[1]]))) x <- ss3_fwd(x)
  ey <- attr(x, "endyr") %||% attr(x[[1]], "endyr")
  start <- start %||% (ey - 10)
  # reference (end) year and implementation (advice) year, as the FLRef summary sheets: the intermediate years lie
  # between the two lines
  vlines <- vlines %||% unique(c(0, .fwd_nint(x, ey) + 1))
  a <- .win_stk(stk %||% x[[1]], start = start, end = ey)
  attr(a, "endyr") <- ey
  ca <- as.numeric(FLCore::iterMedians(FLCore::catch(a)[, as.character(ey)]))
  sc <- lapply(as.list(x), function(s){
    w <- .win_stk(s, start = ey, end = end); attr(w, "endyr") <- ey
    # catch of the end year from the assessment when the scenario has forecast catch only (SSsummarize input)
    if(is.finite(ca) && all(!is.finite(FLCore::catch(w)[, as.character(ey)]) | FLCore::catch(w)[, as.character(ey)] == 0))
      FLCore::catch(w)[, as.character(ey)] <- ca
    w })
  st <- FLCore::FLStocks(c(list(Assessment = a), sc))
  cols <- c(Assessment = "black", stats::setNames(colours %||% ss3runcols(length(sc)), names(sc)))
  p <- ss3plot_trj(st, type = type, probs = probs, colours = cols, ncol = ncol, fwd = TRUE, ci = ci,
                   refpts = "first", refnames = refnames, ...)
  # replace the generic forecast marker by the advice lines
  p$layers <- Filter(function(l) !inherits(l$geom, "GeomVline"), p$layers)
  p$labels$caption <- NULL
  if(length(vlines)) p <- p + ggplot2::geom_vline(xintercept = ey + vlines, linetype = 2, linewidth = 0.3)
  p + ggplot2::labs(colour = NULL)
  })(), text)
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
#' @param unc optional FLStocks with iterations for the probabilities (default \code{x} when it has iterations;
#'   otherwise lognormal with the MLE standard error of SSB)
#' @param blim Blim (default the \code{Blim} reference point of the stocks)
#' @return data.table with Basis, catch, F, SSB(year + 1), SSB change (\%) and P(SSB < Blim)
#' @examples
#' fw <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18_fwd.rds", package = "ss3diags2"))
#' fs <- ss3_fwd(fw, Fref = "Btgt", bpa = 0.5, blim = 0.25)
#' ss3_fwdtable(fs)
#' ss3_fwdtable(fs, year = 2027, blim = 30000)
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
    se <- attr(s, "se.log")
    # Blim of each draw (ensembles: each model's own Blim), else one value
    bli <- if(is.null(blim) && !is.null(unc) && "Blim" %in% dimnames(FLCore::refpts(unc[[n]]))$params &&
              FLCore::dims(FLCore::refpts(unc[[n]]))$iter > 1) as.numeric(FLCore::refpts(unc[[n]])["Blim", ]) else bl
    pb <- if(!is.null(unc) && all(is.finite(bli))) mean(as.numeric(FLCore::ssb(unc[[n]])[, y1]) < bli) else
      if(!is.null(se) && is.finite(bl)) stats::pnorm((log(bl) - log(b1))/as.numeric(se$SSB[, y1])) else NA_real_
    d <- data.table(Basis = n, C = round(med(FLCore::catch(s)[, y]), 1), F = round(med(FLCore::fbar(s)[, y]), 3),
                    SSB = round(b1), SSB_change = round(100*(b1/b0 - 1), 1), P_Blim = round(pb, 3))
    setnames(d, c("C", "F", "SSB", "P_Blim"), c(paste0("C_", y), paste0("F_", y), paste0("SSB_", y1),
                                                 paste0("P(SSB_", y1, "<Blim)")))
    d
  }))
}


#' Ensemble short-term forecast
#'
#' Combines the forecast scenarios of several assessment models into a model ensemble, scenario by scenario: each
#' model contributes draws in proportion to its weight (as \code{\link{ss3_ensemble}()}, and the ss3diags
#' \code{SSdeltaMVLN(weight = )} + FLRef \code{ss2FLStockR()} workflow). Ratios to the targets are computed within each
#' model, and each draw keeps its model's Blim, so \code{\link{ss3_fwdtable}()} gives the ensemble risk
#' P(SSB < Blim) as the share of draws below their own model's Blim.
#'
#' @param x named list of models, each the forecast scenarios of one model: FLStocks from
#'   \code{ss3_fwd(output = "iters")}, a named list of SS3 reports (forecast runs) or their \code{r4ss::SSsummarize()}
#'   summary (lognormal draws from the MLE and SEs; exact for the marginal risk in a year). All models need the same
#'   scenario names
#' @param weights model weights (default equal); normalised to sum to one
#' @param iters average number of draws per model (see \code{\link{ss3_ensemble}()})
#' @param seed random seed
#' @param ... further arguments to \code{\link{ss3_fwd}()} for SS3 reports (e.g. \code{Fref}, \code{blim}, \code{bpa})
#' @return FLStocks of the joint ensemble by scenario, with attributes \code{endyr}, \code{weights} and
#'   \code{models} (the resampled FLStocks of each model by scenario)
#' @seealso \code{\link{ss3plot_fwd}()}, \code{\link{ss3_fwdtable}()}
#' @examples
#' \donttest{
#' # two models (here the same run twice, with different seeds) and one scenario
#' hke <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18.rds", package = "ss3diags2"))
#' fa <- ss3_fwd(list(Fadv = hke), output = "iters", iters = 200, Fref = "Btgt", bpa = 0.5, blim = 0.25)
#' fb <- ss3_fwd(list(Fadv = hke), output = "iters", iters = 200, Fref = "Btgt", bpa = 0.5, blim = 0.2)
#' fe <- ss3_fwdensemble(list(A = fa, B = fb), weights = c(0.6, 0.4), iters = 200)
#' ss3_fwdtable(fe)
#' ss3plot_fwd(fe, start = 2010)
#' }
#' @export
ss3_fwdensemble <- function(x, weights = NULL, iters = 1000, seed = 1, ...){
  x <- lapply(x, function(m) if(.is_summary(m)) ss3_fwd(m, output = "iters", iters = iters, ...) else m)
  if(!is.list(x) || length(x) < 2) stop("x must be a named list of at least two models", call. = FALSE)
  nm <- names(x); if(is.null(nm) || any(nm == "")) nm <- paste0("model", seq_along(x))
  names(x) <- nm
  x <- lapply(x, function(m){
    if(is(m, "FLStocks")) return(m)
    if(is.list(m) && length(m) && all(vapply(m, .is_ss3rep, TRUE))) return(ss3_fwd(m, output = "iters", iters = iters, ...))
    stop("each model must be FLStocks from ss3_fwd(output = \"iters\") or a list of SS3 forecast reports",
         call. = FALSE)
  })
  sc <- Reduce(intersect, lapply(x, names))
  if(!length(sc)) stop("the models have no scenario names in common", call. = FALSE)
  miss <- setdiff(unique(unlist(lapply(x, names))), sc)
  if(length(miss)) warning("scenarios not in all models dropped: ", paste(miss, collapse = ", "), call. = FALSE)
  ens <- lapply(sc, function(s) suppressWarnings(ss3_ensemble(lapply(x, `[[`, s), weights = weights, iters = iters,
                                                              seed = seed)))
  out <- FLCore::FLStocks(stats::setNames(lapply(ens, function(e) e[["Joint"]]), sc))
  for(s in sc) attr(out[[s]], "endyr") <- min(vapply(x, function(m) attr(m[[s]], "endyr") %||% NA_real_, 1), na.rm = TRUE)
  attr(out, "endyr") <- min(vapply(x, function(m) attr(m, "endyr") %||% attr(m[[1]], "endyr") %||% NA_real_, 1),
                            na.rm = TRUE)
  attr(out, "weights") <- attr(ens[[1]], "weights")
  attr(out, "models") <- stats::setNames(lapply(ens, function(e) e[setdiff(names(e), "Joint")]), sc)
  # intermediate years from the models' MLE (the pooled draws carry no MLE)
  attr(out, "nint") <- .fwd_nint(x[[1]][sc], attr(out, "endyr"))
  for(s in sc){ m <- attr(x[[1]][[s]], "mle"); if(!is.null(m)) attr(out[[s]], "mle") <- m }
  out
}


#' Lognormal draws from an MLE stock with standard errors
#'
#' For forecast runs available only as \code{r4ss::SSsummarize()} output (no covariance): draws of SSB, F and
#' recruitment from the MLE and the log-scale SEs (\code{attr(s, "se.log")}). One standard normal deviate per draw and
#' quantity for all years, so the marginal distribution in each year (e.g. P(SSB < Blim)) is exact; correlations
#' between years and quantities are not represented. Catch is kept at the MLE (fixed in catch scenarios).
#'
#' @param s FLStockR (MLE) with \code{attr(s, "se.log")}
#' @param iters number of draws
#' @param seed random seed
#' @return FLStockR with iterations, the attributes of \code{s} and the B/Btgt and F/Ftgt draws in \code{status}
#' @noRd
.lnorm_iters <- function(s, iters, seed = 1){
  se <- attr(s, "se.log")
  if(is.null(se)) stop("no standard errors: draws need the SEs of the derived quantities", call. = FALSE)
  old <- if(exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if(!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed)
  yrs <- as.numeric(dimnames(s)$year)
  v <- function(q) as.numeric(q)[seq_along(yrs)]
  sev <- function(k){ x <- v(se[[k]][, as.character(yrs)]); x[!is.finite(x)] <- 0; x }
  z <- matrix(stats::rnorm(3*iters), ncol = 3)
  dr <- function(m, sd, zz) exp(log(pmax(m, 1e-12)) - sd^2/2*0 + outer(sd, zz))     # median-unbiased (MLE = median)
  ssb <- dr(v(FLCore::ssb(s)), sev("SSB"), z[, 1]); fb <- dr(v(FLCore::fbar(s)), sev("F"), z[, 2])
  rec <- dr(v(FLCore::rec(s)), sev("Rec"), z[, 3])
  d <- data.table(year = rep(yrs, iters), iter = rep(seq_len(iters), each = length(yrs)), SSB = c(ssb), F = c(fb),
                  Recr = c(rec), Catch = rep(v(FLCore::catch(s)), iters))
  out <- FLCore::FLStockR(.stk_build(d, FLCore::name(s)), refpts = FLCore::refpts(s))
  for(a in c("basis", "endyr", "se.log")) attr(out, a) <- attr(s, a)
  # MLE trajectories (used e.g. to find the intermediate years, where the draws would differ by chance)
  bt0 <- .rp(s, "Btgt"); ft0 <- .rp(s, "Ftgt")
  attr(out, "mle") <- data.table(year = yrs, SSB = v(FLCore::ssb(s)), F = v(FLCore::fbar(s)), Recr = v(FLCore::rec(s)),
                                 Catch = v(FLCore::catch(s)), stock = v(FLCore::ssb(s))/bt0,
                                 harvest = v(FLCore::fbar(s))/ft0, type = "fit")
  bt <- .rp(s, "Btgt"); ft <- .rp(s, "Ftgt")
  mk <- function(x) long2flq(data.table(age = "1", year = rep(yrs, iters), iter = rep(seq_len(iters), each = length(yrs)),
                                        data = c(x)), quant = "age", units = "")
  attr(out, "status") <- list(stock = mk(ssb/bt), harvest = mk(fb/ft))
  out
}

#' Risk and catch by catch scenario
#'
#' Fishing opportunities over a grid of catch (TAC) scenarios: catch, F, SSB and P(SSB < Blim) in the advice year by
#' scenario, for each model and the ensemble (from \code{\link{ss3_fwdensemble}()}), and the catch at which the risk
#' reaches \code{risk} (linear interpolation between neighbouring scenarios).
#'
#' @param x FLStocks of scenarios: \code{\link{ss3_fwdensemble}()} (ensemble and models) or \code{ss3_fwd(output =
#'   "iters")} (one model)
#' @param year advice year (catch year; SSB at the start of the next year); default as \code{\link{ss3_fwdtable}()}
#' @param risk target risk P(SSB < Blim)
#' @param value numeric value of each scenario (e.g. the TAC or multiplier); default the number in the scenario
#'   names (TAC_1.25 -> 1.25)
#' @return list: \code{table} (data.table model, scenario, value, catch, F, SSB, P(SSB < Blim)) and \code{at.risk}
#'   (value and catch at the target risk by model)
#' @seealso \code{\link{ss3plot_fwdrisk}()}
#' @examples
#' \donttest{
#' fw <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18_fwd.rds", package = "ss3diags2"))
#' fs <- ss3_fwd(fw, output = "iters", iters = 500, Fref = "Btgt", blim = 0.6)
#' ss3_fwdrisk(fs, value = c(1, 0.85, 0))
#' }
#' @export
ss3_fwdrisk <- function(x, year = NULL, risk = 0.05, value = NULL){
  ey <- attr(x, "endyr") %||% attr(x[[1]], "endyr")
  y <- as.character(year %||% (ey + .fwd_nint(x, ey) + 1)); y1 <- as.character(as.numeric(y) + 1)
  val <- value %||% suppressWarnings(as.numeric(sub("^[^0-9.-]*", "", names(x))))
  if(length(val) != length(x) || any(!is.finite(val))) val <- seq_along(x)
  one <- function(st, model){
    rbindlist(lapply(seq_along(st), function(i){
      s <- st[[i]]
      b <- as.numeric(FLCore::ssb(s)[, y1])
      rp <- FLCore::refpts(s)
      bl <- if("Blim" %in% dimnames(rp)$params) as.numeric(rp["Blim", ]) else NA_real_
      data.table(model = model, scenario = names(st)[i], value = val[i],
                 catch = stats::median(as.numeric(FLCore::catch(s)[, y])),
                 F = stats::median(as.numeric(FLCore::fbar(s)[, y])), SSB = stats::median(b),
                 risk = if(all(is.finite(bl))) mean(b < bl) else NA_real_)
    }))
  }
  mods <- attr(x, "models")
  tb <- one(x, if(is.null(mods)) "model" else "Ensemble")
  if(!is.null(mods)){
    mn <- names(mods[[1]])
    tb <- rbind(tb, rbindlist(lapply(mn, function(m)
      one(FLCore::FLStocks(stats::setNames(lapply(names(x), function(s) mods[[s]][[m]]), names(x))), m))))
  }
  setorder(tb, model, value)
  target <- risk
  at <- rbindlist(lapply(split(tb, by = "model"), function(d){
    d <- d[order(value)]; r <- d$risk
    j <- which(diff(r >= target) != 0)                  # neighbouring scenarios on both sides of the target risk
    if(!length(j) || anyNA(r[c(j[1], j[1] + 1)])) return(data.table(model = d$model[1], value = NA_real_,
                                                                      catch = NA_real_))
    j <- j[1]; f <- (target - r[j])/(r[j + 1] - r[j])
    data.table(model = d$model[1], value = d$value[j] + f*(d$value[j + 1] - d$value[j]),
               catch = d$catch[j] + f*(d$catch[j + 1] - d$catch[j]))
  }))
  list(table = tb[], at.risk = at[], risk = risk, year = as.numeric(y))
}

#' Plot risk by catch scenario
#'
#' P(SSB < Blim) in the year after the advice year against the catch (TAC) scenario, for each model and the ensemble,
#' with the target risk and the catch at that risk.
#'
#' @param x output of \code{\link{ss3_fwdrisk}()}
#' @param xvar x axis: \code{"catch"} (catch in the advice year) or \code{"value"} (scenario value, e.g. TAC multiplier)
#' @param colours model colours (default r4ss colours; the ensemble in black)
#' @param text show the title, subtitle and caption (\code{FALSE}: none; labels inside the plot and legends stay)
#' @return ggplot
#' @examples
#' \donttest{
#' fw <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18_fwd.rds", package = "ss3diags2"))
#' fs <- ss3_fwd(fw, output = "iters", iters = 500, Fref = "Btgt", blim = 0.6)
#' ss3plot_fwdrisk(ss3_fwdrisk(fs, value = c(1, 0.85, 0)))
#' }
#' @export
ss3plot_fwdrisk <- function(x, xvar = c("catch", "value"), colours = NULL, text = TRUE){
  xvar <- match.arg(xvar)
  tb <- copy(x$table); at <- copy(x$at.risk)
  ms <- unique(tb$model); oth <- setdiff(ms, "Ensemble")
  cols <- colours %||% c(stats::setNames(ss3runcols(max(length(oth), 1))[seq_along(oth)], oth), Ensemble = "black")
  tb[, model := factor(model, levels = c(oth, intersect("Ensemble", ms)))]
  at[, model := factor(model, levels = levels(tb$model))]
  tb[, xx := get(xvar)]; at[, xx := get(xvar)]
  p <- ggplot2::ggplot(tb, ggplot2::aes(xx, risk, colour = model)) +
    ggplot2::geom_hline(yintercept = x$risk, linetype = 2, colour = "red") +
    ggplot2::geom_line(ggplot2::aes(linewidth = model == "Ensemble")) + ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_linewidth_manual(values = c(`FALSE` = 0.6, `TRUE` = 1.2), guide = "none") +
    ggplot2::scale_colour_manual(values = cols, name = NULL) +
    ggplot2::labs(x = if(xvar == "catch") paste("Catch", x$year) else "Scenario", y = paste0("P(SSB", x$year + 1, " < Blim)"),
                  title = "Risk by catch scenario",
                  subtitle = paste0("dashed: target risk ", x$risk, "; lines: catch at the target risk")) +
    theme_ss3()
  if(any(is.finite(at$xx))) p <- p + ggplot2::geom_vline(data = at[is.finite(xx)], ggplot2::aes(xintercept = xx, colour = model),
                                                          linetype = 3)
  .drop_text(p, text)
}
