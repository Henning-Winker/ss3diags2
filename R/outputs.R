# outputs.R - model outputs: reference points, MLE and MVLN trajectories, FLStockR, equilibrium curve

# reference point basis of an SS3 run (logic of ss3diags::SSdeltaMVLN)
.ss3_basis <- function(ss3rep, Fref = NULL){
  hat <- as.data.table(ss3rep$derived_quants)
  lab <- .col(hat, c("Label", "LABEL")); val <- .col(hat, c("Value", "VALUE"))
  g <- function(l){ v <- val[lab %in% l]; if(length(v)) v[1] else NA_real_ }
  endyr <- ss3rep$endyr
  bt <- g(paste0("SSB_", endyr)); bv <- g(c("SSB_virgin", "SSB_Virgin")); b0 <- g(c("SSB_unfished", "SSB_Unfished"))
  fmsy <- g(c("Fstd_MSY", "annF_MSY")); fspr <- g(c("Fstd_SPR", "annF_SPR", "Fstd_SPRtgt")); bmsy <- g("SSB_MSY")
  btgt <- g("SSB_Btgt"); f01 <- "SSB_F01" %in% lab; bf01 <- FALSE
  bspr <- g(c("SSB_SPR", "SSB_SPRtgt")); fbtgt <- g(c("Fstd_Btgt", "annF_Btgt"))
  # same benchmark (MSY option set to F(Btgt) or F(SPR)): only for estimated (non-zero) quantities
  near <- function(a, b) isTRUE(is.finite(a) && is.finite(b) && a > 0 && b > 0 && abs(a/b - 1) < 1e-4)
  if(f01){
    if(is.null(Fref)) Fref <- "Btgt"
    btgt <- g("SSB_F01")*b0
    if(round(g("annF_F01"), 2) == round(fmsy, 2)) bf01 <- TRUE
  }
  bratio <- g(paste0("Bratio_", endyr))
  bb.check <- c(bt/bv, bt/bmsy, bt/btgt)
  # as SSdeltaMVLN: MSY benchmark equal to the Btgt or SPR benchmark sets the F reference
  option.btgt <- near(btgt, bmsy)
  if(option.btgt && is.null(Fref)) Fref <- "Btgt"
  if(near(fmsy, fspr) && is.null(Fref)) Fref <- "SPR"
  bb <- max(which(abs(bratio - bb.check) == min(abs(bratio - bb.check), na.rm = TRUE)))
  if(bf01) bb <- 4
  fbasis <- ss3rep$F_std_basis %||% ss3rep$F_report_basis
  if(is.null(fbasis)) stop("F_report_basis not found in ss3rep", call. = FALSE)
  fbasis <- strsplit(fbasis, ";")[[1]][1]
  btarg <- if(is.na(ss3rep$btarg %||% NA)) 0 else ss3rep$btarg
  # F basis (pattern match to cover SS3 versions): 1 absolute F, 2 F/Fmsy, 3 F/F at Btgt, 4 F/F at SPR target
  fb <- if(grepl("abs_F", fbasis)) 1L else if(grepl("SPR", fbasis)) 4L else if(grepl("F_at_B|B[0-9]+%", fbasis)) 3L else
    if(grepl("msy", fbasis, ignore.case = TRUE)) 2L else NA_integer_
  if(is.na(fb)) stop("F_report_basis '", fbasis, "' not supported; rerun SS3 with F_report_basis = 1 in starter.ss",
                     call. = FALSE)
  if(is.null(Fref)) Fref <- c("MSY", "MSY", "Btgt", "SPR")[fb]
  if(Fref == "F01") Fref <- "Btgt"
  ok <- (fb == 1) || (fb == 2 && Fref == "MSY") || (fb == 3 && Fref == "Btgt") || (fb == 4 && Fref == "SPR")
  if(!ok) stop("Fref = ", Fref, " conflicts with F_report_basis ", fbasis, call. = FALSE)
  Fquant <- Fref; if(Fquant == "Btgt" && f01) Fquant <- "F01"
  Fname <- if("Fstd_MSY" %in% lab) "Fstd_" else "annF_"
  bref <- if(btarg > 0) btarg else round(btgt/b0, 2)
  if(!is.finite(bref) || bref <= 0) bref <- 0.4                           # as SSsettingsBratioF
  trg <- round(bref*100)
  sprt <- ss3rep$sprtarg %||% NA; if(!is.finite(sprt) || sprt <= 0) sprt <- 0.4
  spr <- round(sprt*100)
  # labels: name the benchmark behind the target (MSY option may be set to F(Btgt) or F(SPR))
  blab <- switch(as.character(bb), "1" = paste0("B", trg), "3" = paste0("B", trg), "4" = "BF0.1",
                 "2" = if(near(bmsy, btgt)) paste0("B", trg) else if(near(bmsy, bspr)) paste0("Bspr", spr) else "Bmsy")
  flab <- switch(Fquant, Btgt = paste0("F", trg), SPR = paste0("Fspr", spr), F01 = "F0.1",
                 MSY = if(near(fmsy, fbtgt)) paste0("F", trg) else if(near(fmsy, fspr)) paste0("Fspr", spr) else "Fmsy")
  proxy <- near(bmsy, btgt) || near(bmsy, bspr) || near(fmsy, fspr) || near(fmsy, fbtgt)
  clab <- if(!proxy) "MSY" else if(near(bmsy, bspr) || near(fmsy, fspr)) paste0("Yspr", spr) else paste0("Y", trg)
  list(bb = bb, fb = fb, Fquant = Fquant, Fname = Fname, bref = bref, bv = bv, b0 = b0, option.btgt = option.btgt,
       bbasis = c("SSB/SSB0", "SSB/SSBMSY", "SSB/SSBtgt", "SSB/SSBF01")[bb], fbasis = fbasis,
       Bref = blab, Fref = flab, Cref = clab, proxy = proxy,
       msy = g(c("Dead_Catch_MSY", "TotYield_MSY", "Ret_Catch_MSY")), r0 = g(c("Recr_unfished", "Recr_Unfished", "Recr_Virgin")))
}

# MLE trajectories from the time series when derived quantities are not available (run without hessian/sdreport)
.ss3_ts_mle <- function(ss3rep, btgt = NULL){
  ts <- as.data.table(ss3rep$timeseries)
  tt <- ts[Era == "TIME"]
  ssb <- tt[, .(SSB = sum(SpawnBio, na.rm = TRUE)), by = .(year = Yr)]
  rec <- tt[, .(Recr = sum(Recruit_0, na.rm = TRUE)), by = .(year = Yr)]
  ct <- as.data.table(ss3rep$catch)
  dcol <- intersect(c("dead_bio", "kill_bio", "Exp"), names(ct))[1]
  cat <- ct[Yr %in% tt$Yr, .(Catch = sum(get(dcol), na.rm = TRUE)), by = .(year = Yr)]
  ex <- as.data.table(ss3rep$exploitation)
  fn <- intersect(.fleetnames(ss3rep), names(ex))
  Fy <- if(nrow(ex) && "annual_F" %in% names(ex) && any(ex$annual_F > 0, na.rm = TRUE)) {
    ex[, .(F = sum(annual_F*Seas_dur, na.rm = TRUE)), by = .(year = Yr)]
  } else if(nrow(ex) && length(fn)){
    message("annual F not reported (run without sdreport?): F = sum of fleet apical F")
    ex[, .(F = sum(rowSums(.SD, na.rm = TRUE)*Seas_dur)), by = .(year = Yr), .SDcols = fn]
  } else data.table(year = ssb$year, F = NA_real_)
  m <- Reduce(function(a, b) merge(a, b, by = "year", all.x = TRUE), list(ssb, rec, cat, Fy))
  v <- ts[Era == "VIRG"]
  b0 <- sum(v$SpawnBio, na.rm = TRUE); r0 <- sum(v$Recruit_0, na.rm = TRUE)
  bt <- btgt %||% (if(isTRUE(ss3rep$btarg > 0)) ss3rep$btarg else NA_real_)
  m[, `:=`(type = "fit", stock = SSB/(bt*b0), harvest = NA_real_)]
  list(mle = m[year <= ss3rep$endyr], refpts = FLCore::FLPar(Ftgt = NA_real_, Btgt = bt*b0, MSY = NA_real_, B0 = b0, R0 = r0),
       basis = list(B = if(is.finite(bt)) paste0("B", round(100*bt)) else "Btgt", F = "Ftgt", bbasis = "timeseries",
                    fbasis = "timeseries", Btgtref = bt, ssbunits = ss3rep$SpawnOutputUnits %||% "biomass"))
}

.derived_ok <- function(ss3rep){
  h <- ss3rep$derived_quants; if(is.null(h)) return(FALSE)
  lab <- .col(h, c("Label", "LABEL")); val <- .col(h, c("Value", "VALUE"))
  v <- val[lab == paste0("SSB_", ss3rep$endyr)]
  length(v) == 1 && isTRUE(v > 0)
}

#' Equilibrium yield curve
#'
#' The SS3 equilibrium yield curve (\code{ss3rep$equil_yield}, the SPR/YPR profile): SSB, F and yield across the range
#' of F explored by SS3, with the MSY reference points, for production function plots and \code{\link{ss3_spm}()}.
#'
#' The SS3 table stacks several profiles (\code{SPRloop}): repeated F grids from F = 0 (loops 0 to 2), single rows
#' for the benchmarks (SPR target, B target, MSY, ...) and, in some versions, an extension from Fmsy to high F
#' (loop 7). By default the last F grid from F = 0 (the converged loop, which reproduces the SS3 MSY) is combined with
#' the extension loops; earlier grids can differ (e.g. another F scale in loop 0).
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param loops optional SPR loops to use (default: the last grid starting at F = 0 plus extension loops)
#' @return data.table with columns \code{F}, \code{SSB}, \code{yield}, \code{loop} (sorted by SSB); attributes
#'   \code{Fmsy}, \code{Bmsy}, \code{MSY}, \code{B0}, \code{sd} (SS3 standard errors of MSY and B0) and
#'   \code{ssbunits}
#' @export
ss3_eqcurve <- function(ss3rep, loops = NULL){
  eq <- as.data.table(ss3rep$equil_yield)
  if(!nrow(eq)) stop("no equil_yield table in ss3rep", call. = FALSE)
  d <- data.table(loop = .col(eq, "SPRloop", 0), F = .col(eq, c("F_report", "F_std")), SSB = eq$SSB,
                  yield = .col(eq, c("Tot_Catch", "Catch")))
  d <- d[is.finite(F) & is.finite(SSB) & is.finite(yield)]
  if(is.null(loops)){
    ln <- d[, .(n = .N, f0 = min(F), f1 = max(F)), by = loop][n > 2]
    grid <- ln[f0 <= 0.01*f1]$loop; ext <- ln[f0 > 0.01*f1]$loop
    loops <- c(if(length(grid)) max(grid), ext)
  }
  d <- d[loop %in% loops][order(SSB, -F)]
  d <- unique(d, by = c("F", "SSB"))
  hat <- as.data.table(ss3rep$derived_quants)
  lab <- .col(hat, c("Label", "LABEL")); val <- .col(hat, c("Value", "VALUE")); sdv <- .col(hat, c("StdDev", "STDDEV"))
  g <- function(l, v = val){ x <- v[lab %in% l]; if(length(x)) x[1] else NA_real_ }
  ms <- c("Dead_Catch_MSY", "Ret_Catch_MSY"); bs <- c("SSB_unfished", "SSB_Unfished")
  structure(d, Fmsy = g(c("Fstd_MSY", "annF_MSY")), Bmsy = g("SSB_MSY"), MSY = g(ms), B0 = g(bs),
            sd = c(MSY = g(ms, sdv), B0 = g(bs, sdv)), ssbunits = ss3rep$SpawnOutputUnits %||% "biomass")
}

utils::globalVariables(c("label.i", "label.j", "corr", "stock", "harvest", "SSB", "Recr", "Catch", "F", "Era",
                         "SpawnBio", "Recruit_0", "annual_F", "Seas_dur"))
