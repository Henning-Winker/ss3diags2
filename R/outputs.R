# outputs.R - model outputs: reference points, MLE and MVLN trajectories, FLStockR, equilibrium curve

# reference point basis of an SS3 run (logic of ss3diags::SSdeltaMVLN)
.ss3_basis <- function(ss3rep, Fref = NULL){
  hat <- as.data.table(ss3rep$derived_quants)
  lab <- .col(hat, c("Label", "LABEL")); val <- .col(hat, c("Value", "VALUE"))
  g <- function(l){ v <- val[lab %in% l]; if(length(v)) v[1] else NA_real_ }
  endyr <- ss3rep$endyr
  bt <- g(paste0("SSB_", endyr)); bv <- g(c("SSB_virgin", "SSB_Virgin")); b0 <- g(c("SSB_unfished", "SSB_Unfished"))
  fmsy <- g(c("Fstd_MSY", "annF_MSY")); fspr <- g(c("Fstd_SPR", "annF_SPR")); bmsy <- g("SSB_MSY")
  btgt <- g("SSB_Btgt"); f01 <- "SSB_F01" %in% lab; bf01 <- FALSE
  if(f01){
    if(is.null(Fref)) Fref <- "Btgt"
    btgt <- g("SSB_F01")*b0
    if(round(g("annF_F01"), 2) == round(fmsy, 2)) bf01 <- TRUE
  }
  bratio <- g(paste0("Bratio_", endyr))
  bb.check <- c(bt/bv, bt/bmsy, bt/btgt)
  option.btgt <- isTRUE(all.equal(btgt, bmsy))
  if(option.btgt && is.null(Fref)) Fref <- "Btgt"
  if(isTRUE(all.equal(fmsy, fspr)) && is.null(Fref)) Fref <- "SPR"
  bb <- max(which(abs(bratio - bb.check) == min(abs(bratio - bb.check), na.rm = TRUE)))
  if(bf01) bb <- 4
  fbasis <- ss3rep$F_std_basis %||% ss3rep$F_report_basis
  if(is.null(fbasis)) stop("F_report_basis not found in ss3rep", call. = FALSE)
  fbasis <- strsplit(fbasis, ";")[[1]][1]
  btarg <- if(is.na(ss3rep$btarg %||% NA)) 0 else ss3rep$btarg
  fopts <- c("_abs_F", "(F)/(Fmsy)", paste0("(F)/(F_at_B", btarg*100, "%)"), paste0("(F)/(F", btarg*100, "%SPR)"))
  fb <- match(fbasis, fopts)
  if(is.na(fb)) stop("F_report_basis '", fbasis, "' not supported; rerun SS3 with F_report_basis = 1 in starter.ss",
                     call. = FALSE)
  if(is.null(Fref)) Fref <- c("MSY", "MSY", "Btgt", "SPR")[fb]
  if(Fref == "F01") Fref <- "Btgt"
  ok <- (fb == 1) || (fb == 2 && Fref == "MSY") || (fb == 3 && Fref == "Btgt") || (fb == 4 && Fref == "SPR")
  if(!ok) stop("Fref = ", Fref, " conflicts with F_report_basis ", fbasis, call. = FALSE)
  Fquant <- Fref; if(Fquant == "Btgt" && f01) Fquant <- "F01"
  Fname <- if("Fstd_MSY" %in% lab) "Fstd_" else "annF_"
  bref <- if(btarg > 0) btarg else round(btgt/b0, 2)
  trg <- round(bref*100); spr <- round((ss3rep$sprtarg %||% 0.4)*100)
  list(bb = bb, fb = fb, Fquant = Fquant, Fname = Fname, bref = bref, bv = bv, b0 = b0, option.btgt = option.btgt,
       bbasis = c("SSB/SSB0", "SSB/SSBMSY", "SSB/SSBtgt", "SSB/SSBF01")[bb], fbasis = fbasis,
       Bref = c(paste0("B", trg), "Bmsy", paste0("B", trg), "BF0.1")[bb],
       Fref = c(MSY = "Fmsy", Btgt = paste0("Fb", trg), SPR = paste0("F", spr), F01 = "F0.1")[[Fquant]],
       msy = g(c("Dead_Catch_MSY", "Ret_Catch_MSY")), r0 = g(c("Recr_unfished", "Recr_Unfished", "Recr_Virgin")))
}

#' Multivariate lognormal (MVLN) approximation of SS3 stock status
#'
#' Draws from the joint distribution of B/Bref, F/Fref and Fref for each year from the SS3 estimates and their
#' correlations (CoVar), with SSB and recruitment drawn from their marginal lognormals, as
#' \code{ss3diags::SSdeltaMVLN()} (Walter and Winker 2019), but without dplyr and returning long tables.
#'
#' @param ss3rep output of \code{r4ss::SS_output()} (with \code{covar = TRUE})
#' @param Fref reference for F: \code{"MSY"}, \code{"Btgt"}, \code{"SPR"} or \code{"F01"} (default from the SS3
#'   F_report_basis)
#' @param years years (default: years up to the end year; \code{forecast = TRUE} adds forecast years)
#' @param mc number of draws
#' @param forecast include forecast years
#' @param virgin B/B0 relative to virgin (rather than unfished) SSB
#' @param catch catch column of \code{ss3rep$catch} summed over fleets (\code{"Exp"} default, \code{"Obs"},
#'   \code{"kill_bio"}, \code{"ret_bio"}); forecast catch from \code{ForeCatch_}
#' @param seed random seed (set for each year, as SSdeltaMVLN)
#' @return list of class \code{ss3mvln}: \code{kb} (data.table: year, type, iter, stock = B/Bref, harvest = F/Fref,
#'   SSB, F, Recr, Catch), \code{mle} (the same at the estimates), \code{refpts} (FLPar: Ftgt, Btgt, MSY, B0, R0),
#'   \code{basis} (B and F reference labels)
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' mv <- ss3_mvln(ss3sma, mc = 1000)
#' mv$refpts
#' }
#' @export
ss3_mvln <- function(ss3rep, Fref = NULL, years = NULL, mc = 1000, forecast = FALSE, virgin = FALSE,
                     catch = c("Exp", "Obs", "kill_bio", "ret_bio"), seed = 123){
  catch <- match.arg(catch)
  .need("mvtnorm")
  if(is.null(ss3rep$CoVar)) stop("CoVar required: read the SS3 run with r4ss::SS_output(covar = TRUE)", call. = FALSE)
  b <- .ss3_basis(ss3rep, Fref)
  hat <- as.data.table(ss3rep$derived_quants)
  lab <- .col(hat, c("Label", "LABEL")); val <- .col(hat, c("Value", "VALUE")); sdv <- .col(hat, c("StdDev", "STDDEV"))
  hv <- function(l) val[match(l, lab)]; hs <- function(l) sdv[match(l, lab)]
  allyrs <- suppressWarnings(as.numeric(sub("Bratio_", "", lab[grep("^Bratio_", lab)])))
  allyrs <- allyrs[!is.na(allyrs)][-1]
  yrs <- if(!is.null(years)) years[years %in% allyrs] else if(forecast) allyrs else allyrs[allyrs <= ss3rep$endyr]
  cv <- as.data.table(ss3rep$CoVar)
  cv[label.j == "_", label.j := label.i]
  fr <- paste0(b$Fname, b$Fquant)
  fref <- c(hv(fr), hs(fr))
  sy <- paste0(rep(c("Bratio_", "F_"), each = length(yrs)), yrs)       # correlations within the selected years only
  cvs <- cv[label.i %in% sy | (label.i == fr & label.j %in% sy)]
  mcor <- function(i, j) { x <- cvs[label.i %in% i & label.j %in% j, corr]; if(length(x)) mean(x) else 0 }
  out <- lapply(yrs, function(yr){
    set.seed(seed)
    ny <- c(yr - 1, yr, yr + 1)
    fv <- hv(paste0("F_", yr)); fv <- ifelse(fv == 0, 0.001, fv); fs <- hs(paste0("F_", yr))
    bv <- hv(paste0("Bratio_", yr)); bs <- hs(paste0("Bratio_", yr))
    varF <- log(1 + (fs/fv)^2); varB <- log(1 + (bs/bv)^2); varR <- log(1 + (fref[2]/fref[1])^2)
    cBF <- log(1 + mcor(paste0("Bratio_", ny), paste0("F_", ny))*sqrt(varF*varB))
    cBR <- log(1 + mcor(fr, paste0("Bratio_", ny))*sqrt(varB*varR))
    cFR <- log(1 + mcor(fr, paste0("F_", ny))*sqrt(varF*varR))
    S <- matrix(c(varB, cBF, cFR, cBF, varF, cBR, cFR, cBR, varR), 3)   # as SSdeltaMVLN: [2,3] = cov(B, Fref)
    x <- exp(mvtnorm::rmvnorm(mc, mean = log(c(bv, fv, fref[1])), sigma = S, method = "svd"))
    q <- vapply(c("SSB", "Recr"), function(k){
      v <- hv(paste0(k, "_", yr)); s <- sqrt(log(1 + (hs(paste0(k, "_", yr))/v)^2))
      stats::rlnorm(mc, log(v) - 0.5*s*s, s)
    }, numeric(mc))
    list(kb = data.table(year = yr, type = ifelse(yr <= ss3rep$endyr, "fit", "forecast"), iter = seq_len(mc),
                         stock = x[, 1], harvest = x[, 2], F = x[, 3], SSB = q[, "SSB"], Recr = q[, "Recr"]),
         mle = data.table(year = yr, type = ifelse(yr <= ss3rep$endyr, "fit", "forecast"), stock = bv, harvest = fv,
                          F = fref[1], SSB = hv(paste0("SSB_", yr)), Recr = hv(paste0("Recr_", yr))))
  })
  kb <- rbindlist(lapply(out, `[[`, "kb")); mle <- rbindlist(lapply(out, `[[`, "mle"))
  # rescale B ratio to the reference, F absolute and relative, catch (as SSdeltaMVLN)
  rs <- 1
  if(!virgin && (b$bb == 1 || (b$bb == 3 && !b$option.btgt))) rs <- b$bv/b$b0
  if(virgin && b$bb == 3 && b$option.btgt) rs <- b$b0/b$bv
  if(b$bb == 1) rs <- rs/b$bref
  ct <- as.data.table(ss3rep$catch)
  C <- ct[, .(C = sum(get(catch), na.rm = TRUE)), by = .(year = Yr)]
  fy <- suppressWarnings(as.numeric(sub("ForeCatch_", "", lab[grep("^ForeCatch_", lab)])))
  if(length(fy)) C <- rbind(C[!year %in% fy], data.table(year = fy, C = hv(paste0("ForeCatch_", fy))))
  fix <- function(d){
    d[, stock := stock*rs]
    if(b$fb > 1) d[, F := F*harvest] else d[, c("harvest", "F") := list(harvest/F, harvest)]
    d[, Catch := C$C[match(year, C$year)]]
    d
  }
  kb <- fix(kb); mle <- fix(mle)
  ftgt <- mle$F[1]/mle$harvest[1]; btgt <- mle$SSB[1]/mle$stock[1]
  refpts <- FLCore::FLPar(Ftgt = ftgt, Btgt = btgt, MSY = b$msy, B0 = b$b0, R0 = b$r0)
  structure(list(kb = kb[], mle = mle[], refpts = refpts,
                 basis = list(B = b$Bref, F = b$Fref, bbasis = b$bbasis, fbasis = b$fbasis, Btgtref = b$bref,
                              ssbunits = ss3rep$SpawnOutputUnits %||% "biomass")),
            class = "ss3mvln")
}

#' SS3 trajectories as FLStockR
#'
#' Builds the simplified (biomass-dynamic, age "1") \code{FLStockR} of \code{FLRef::ss2FLStockR()} from an SS3 report:
#' \code{ssb()} = SSB, \code{fbar()} = F, \code{rec()} = recruitment, \code{catch()} = catch, with reference points
#' (Ftgt, Btgt, MSY, B0, R0, and Blim, Bpa when given) in \code{@refpts}, so the SS3 run plots and compares with
#' JABBA3, FLicc and FLRef outputs. \code{iters = 0} gives the maximum likelihood estimates, with the log-scale SEs
#' of SSB, F, recruitment, B/Btgt and F/Ftgt attached as \code{attr(x, "se.log")} (FLQuants) for Wald-type intervals;
#' \code{iters > 0} gives MVLN draws (\code{\link{ss3_mvln}()}) as iterations, with the MLE in \code{attr(x, "mle")}
#' and the jointly sampled B/Btgt and F/Ftgt in \code{attr(x, "status")} (FLQuants), used by the Kobe plot.
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param iters 0 (MLE, default) or the number of MVLN draws
#' @param Fref,forecast,catch passed to \code{\link{ss3_mvln}()}
#' @param blim,bpa optional Blim and Bpa as fractions of Btgt, added to the reference points
#' @param name stock name
#' @return FLStockR
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' stk <- ss3_FLStockR(ss3sma)               # MLE
#' stki <- ss3_FLStockR(ss3sma, iters = 500) # MVLN
#' refpts(stk)
#' }
#' @export
ss3_FLStockR <- function(ss3rep, iters = 0, Fref = NULL, forecast = FALSE, catch = c("Exp", "Obs", "kill_bio", "ret_bio"),
                         blim = NULL, bpa = NULL, name = "ss3"){
  catch <- match.arg(catch)
  mv <- ss3_mvln(ss3rep, Fref = Fref, mc = max(iters, 2), forecast = forecast, catch = catch)
  d <- if(iters > 0) mv$kb else mv$mle[, iter := 1L]
  mk <- function(v, u) long2flq(d[, .(age = "1", year, iter, data = get(v))], quant = "age", units = u)
  R <- mk("Recr", "1000"); C <- mk("Catch", "t"); H <- mk("F", "f")
  mat <- long2flq(d[, .(age = "1", year, iter, data = SSB/Recr)], quant = "age", units = "")
  one <- FLCore::FLQuant(1, dimnames = dimnames(R)); zero <- one*0
  stk <- FLCore::FLStock(name = name, stock.n = R, catch.n = C, landings.n = C, discards.n = zero, stock.wt = one,
                         catch.wt = one, landings.wt = one, discards.wt = one, mat = mat, m = one*0.0001, harvest = H,
                         m.spwn = zero, harvest.spwn = zero)
  FLCore::units(FLCore::harvest(stk)) <- "f"
  FLCore::catch(stk) <- FLCore::computeCatch(stk); FLCore::landings(stk) <- FLCore::computeLandings(stk)
  FLCore::stock(stk) <- FLCore::computeStock(stk)
  rp <- mv$refpts
  extra <- c(Blim = blim %||% NA, Bpa = bpa %||% NA); extra <- extra[is.finite(extra)]
  if(length(extra)) rp <- FLCore::FLPar(c(as.numeric(rp), extra*as.numeric(rp["Btgt"])),
                                        params = c(dimnames(rp)$params, names(extra)))
  out <- FLCore::FLStockR(stk, refpts = rp)
  attr(out, "basis") <- mv$basis
  if(iters > 0){
    attr(out, "mle") <- mv$mle
    # jointly sampled B/Btgt and F/Ftgt (SSB and F draws are marginal, so their ratios lose the correlation)
    attr(out, "status") <- FLCore::FLQuants(stock = mk("stock", ""), harvest = mk("harvest", ""))
  } else {
    hat <- as.data.table(ss3rep$derived_quants)
    lab <- .col(hat, c("Label", "LABEL")); val <- .col(hat, c("Value", "VALUE")); sdv <- .col(hat, c("StdDev", "STDDEV"))
    yrs <- mv$mle$year
    sel <- function(k){ i <- match(paste0(k, "_", yrs), lab); sqrt(log(1 + (sdv[i]/val[i])^2)) }
    sq <- function(v) long2flq(data.table(age = "1", year = yrs, data = v), quant = "age", units = "")
    attr(out, "se.log") <- FLCore::FLQuants(SSB = sq(sel("SSB")), F = sq(sel("F")), Rec = sq(sel("Recr")),
                                            stock = sq(sel("Bratio")), harvest = sq(sel("F")))
  }
  attr(out, "eqcurve") <- tryCatch(ss3_eqcurve(ss3rep), error = function(e) NULL)
  out
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

utils::globalVariables(c("label.i", "label.j", "corr", "stock", "harvest", "SSB", "Recr", "Catch", "F"))
