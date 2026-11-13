# mvln.R - multivariate lognormal (MVLN) draws of SS3 stock status, FLStockR (ss3stk) and FLQuants (ss3_flqs)

# catch by year: summed fleet column of ss3rep$catch, forecast catch from ForeCatch_
.ss3_catch_year <- function(ss3rep, catch = "Exp"){
  ct <- as.data.table(ss3rep$catch)
  col <- if(catch %in% names(ct)) catch else sub("^kill_", "dead_", catch)
  C <- ct[, .(C = sum(get(col), na.rm = TRUE)), by = .(year = Yr)]
  hat <- ss3rep$derived_quants; lab <- .col(hat, c("Label", "LABEL")); val <- .col(hat, c("Value", "VALUE"))
  fy <- suppressWarnings(as.numeric(sub("ForeCatch_", "", lab[grep("^ForeCatch_", lab)])))
  if(length(fy)) C <- rbind(C[!year %in% fy], data.table(year = fy, C = val[match(paste0("ForeCatch_", fy), lab)]))
  C
}

# labels of the reference points matching the B ratio and F bases (as SSdeltaMVLN)
.ref_labels <- function(b, hatlab){
  pick <- function(x) x[x %in% hatlab][1]
  B0 <- pick(c("SSB_unfished", "SSB_Unfished"))
  Bl <- switch(as.character(b$bb), "2" = "SSB_MSY", "3" = "SSB_Btgt", NA_character_)
  list(B0 = B0, Btgt = if(is.na(Bl) || !Bl %in% hatlab) NA_character_ else Bl, Bmult = if(is.na(Bl)) b$bref else 1,
       Ftgt = paste0(b$Fname, b$Fquant), MSY = pick(c("Dead_Catch_MSY", "TotYield_MSY", "Ret_Catch_MSY")),
       R0 = pick(c("Recr_unfished", "Recr_Unfished", "Recr_Virgin")))
}

# joint MVLN over all years and reference points (ss3om::ssmvln approach, L.T. Kell); mc = 0: MLE only
.mvln_joint <- function(ss3rep, Fref = NULL, years = NULL, mc = 1000, forecast = FALSE, catch = "Exp", seed = 123){
  b <- .ss3_basis(ss3rep, Fref)
  hat <- as.data.table(ss3rep$derived_quants)
  lab <- .col(hat, c("Label", "LABEL")); val <- .col(hat, c("Value", "VALUE")); sdv <- .col(hat, c("StdDev", "STDDEV"))
  allyrs <- suppressWarnings(as.numeric(sub("Bratio_", "", lab[grep("^Bratio_", lab)])))
  allyrs <- allyrs[!is.na(allyrs)]
  yrs <- if(!is.null(years)) years[years %in% allyrs] else if(forecast) allyrs else allyrs[allyrs <= ss3rep$endyr]
  rl <- .ref_labels(b, lab)
  refl <- unique(stats::na.omit(c(rl$Btgt, rl$B0, rl$Ftgt, rl$MSY, rl$R0)))
  labs <- c(paste0("SSB_", yrs), paste0("Recr_", yrs), paste0("F_", yrs), refl)
  if(any(!labs %in% lab)) stop("not in derived_quants: ", paste(setdiff(labs, lab), collapse = ", "), call. = FALSE)
  v <- val[match(labs, lab)]; s <- sdv[match(labs, lab)]
  zero <- !(is.finite(v) & v > 0)
  v0 <- ifelse(zero, 1, v)
  sl <- ifelse(zero | !is.finite(s), 0, sqrt(log(1 + (s/v0)^2)))     # log-scale SD
  X <- if(mc > 0){
    .need("mvtnorm")
    if(is.null(ss3rep$CoVar)) stop("CoVar required for MVLN draws: read SS3 with r4ss::SS_output(covar = TRUE)", call. = FALSE)
    cv <- as.data.table(ss3rep$CoVar)[label.j != "_" & label.i %in% labs & label.j %in% labs]
    R <- diag(length(labs)); dimnames(R) <- list(labs, labs)
    i <- cbind(match(cv$label.i, labs), match(cv$label.j, labs))
    R[i] <- cv$corr; R[i[, 2:1, drop = FALSE]] <- cv$corr                   # CoVar holds one triangle
    S <- R*outer(sl, sl)
    e <- eigen(S, symmetric = TRUE); ev <- pmax(e$values, 1e-12*max(e$values))
    S <- e$vectors %*% (ev*t(e$vectors))                                   # nearest positive definite
    set.seed(seed)
    x <- exp(mvtnorm::rmvnorm(mc, mean = log(v0), sigma = S, method = "eigen"))
    x[, zero] <- 0; colnames(x) <- labs; x
  } else NULL
  ny <- length(yrs)
  build <- function(m, it){
    g <- function(l) m[, l, drop = TRUE]
    btgt <- if(!is.na(rl$Btgt)) g(rl$Btgt) else rl$Bmult*g(rl$B0)
    ftgt <- g(rl$Ftgt)
    SSB <- m[, paste0("SSB_", yrs), drop = FALSE]; Fy <- m[, paste0("F_", yrs), drop = FALSE]
    Fab <- if(b$fb > 1) Fy*ftgt else Fy; H <- if(b$fb > 1) Fy else Fy/ftgt
    kb <- data.table(year = rep(yrs, each = nrow(m)), type = rep(ifelse(yrs <= ss3rep$endyr, "fit", "forecast"), each = nrow(m)),
                     iter = rep(it, ny), stock = as.vector(SSB/btgt), harvest = as.vector(H), F = as.vector(Fab),
                     SSB = as.vector(SSB), Recr = as.vector(m[, paste0("Recr_", yrs), drop = FALSE]))
    rp <- data.table(iter = it, Ftgt = ftgt, Btgt = btgt, MSY = if(!is.na(rl$MSY)) g(rl$MSY) else NA_real_,
                     B0 = g(rl$B0), R0 = if(!is.na(rl$R0)) g(rl$R0) else NA_real_)
    list(kb = kb, rp = rp)
  }
  mm <- matrix(ifelse(zero, 0, v), nrow = 1, dimnames = list(NULL, labs))
  mle <- build(mm, 1L)
  C <- .ss3_catch_year(ss3rep, catch)
  mle$kb[, Catch := C$C[match(year, C$year)]][, iter := NULL]
  out <- list(mle = mle$kb[], rp.mle = mle$rp, kb = NULL, rp = NULL, b = b)
  if(!is.null(X)){
    dr <- build(X, seq_len(mc)); dr$kb[, Catch := C$C[match(year, C$year)]]
    out$kb <- dr$kb[]; out$rp <- dr$rp[]
  }
  out
}

#' Multivariate lognormal (MVLN) draws of SS3 stock status
#'
#' Draws SSB, recruitment, F, B/Btgt and F/Ftgt, and the reference points, from a multivariate lognormal
#' approximation of the SS3 estimates and their correlations (\code{CoVar}), with data.table (no dplyr).
#'
#' \itemize{
#'   \item \code{method = "joint"} (default), the approach of \code{ss3om::ssmvln()} (L.T. Kell): one MVLN over all
#'     years and the reference points (SSB, recruitment and F of every year, Btgt, Ftgt, MSY, B0, R0) with the full
#'     correlation matrix from CoVar, so each iteration is a coherent trajectory and the reference points have
#'     iterations; B/Btgt = SSB/Btgt and F/Ftgt = F/Ftgt within each iteration. Medians at the MLE (log-scale means
#'     at the log MLE), log-scale SDs \eqn{\sqrt{\log(1 + CV^2)}}; the covariance is made positive definite if needed.
#'   \item \code{method = "annual"}: \code{ss3diags::SSdeltaMVLN()} (Walter and Winker 2019): for each year a
#'     trivariate lognormal of B ratio, F and Fref (neighbouring-year correlations), SSB and recruitment drawn
#'     from their marginals; draws are not linked across years.
#' }
#' The reference points follow the SS3 bases (depletion basis and F_report_basis), as in \code{SSdeltaMVLN()}: Btgt is
#' SSB_MSY (basis SSB/SSBMSY), SSB_Btgt (SSB/SSBtgt) or X x SSB_unfished (SSB/SSB0); Ftgt is Fstd_/annF_ MSY, Btgt or
#' SPR. MLE values and reference points are the same for both methods.
#'
#' @param ss3rep output of \code{r4ss::SS_output()} (with \code{covar = TRUE})
#' @param mc number of draws
#' @param method \code{"joint"} (default) or \code{"annual"}
#' @param Fref reference for F: \code{"MSY"}, \code{"Btgt"}, \code{"SPR"} or \code{"F01"} (default from F_report_basis)
#' @param years years (default up to the end year; \code{fwd = TRUE} adds the forecast years)
#' @param fwd include the SS3 forecast years (e.g. SSB at the start of the year after the end year, needed for the
#'   lagged Kobe plot)
#' @param catch catch column of \code{ss3rep$catch} summed over fleets (\code{"Exp"} default, \code{"Obs"},
#'   \code{"kill_bio"}/\code{"dead_bio"}, \code{"ret_bio"}); forecast catch from \code{ForeCatch_}
#' @param seed random seed
#' @param virgin (\code{"annual"} only) B/B0 relative to virgin rather than unfished SSB
#' @return list of class \code{ss3mvln}: \code{kb} (data.table: year, type, iter, stock = B/Btgt, harvest = F/Ftgt, F,
#'   SSB, Recr, Catch), \code{mle} (the same at the MLE), \code{refpts} (FLPar at the MLE: Ftgt, Btgt, MSY, B0, R0),
#'   \code{refpts.iter} (FLPar with iterations, \code{"joint"}), \code{basis} (B and F reference labels), \code{method}
#' @references
#' Kell, L., Rice, J., Courtney, D. and Winker, H. (2023). Approximation of Kobe posteriors from Stock Synthesis for
#' North Atlantic blue shark. Collect. Vol. Sci. Pap. ICCAT, 80(4): 837-858 (SCRS/2023/143).
#'
#' Walter, J. and Winker, H. (2020). Projections to create Kobe 2 Strategy Matrix using the multivariate log-normal
#' approximation for Atlantic yellowfin tuna. ICCAT SCRS/2019/145.
#'
#' Winker, H., Kell, L.T., Fu, D., Sharma, R., Courtney, D., Carvalho, F., Schirripa, M. and Walter, J. (2019). A rapid
#' approach to approximate Kobe posteriors from age-structured assessment models with applications to North Atlantic
#' shortfin mako. ICCAT SCRS/2019/093.
#' @author Laurence T. Kell (joint MVLN, ss3om), Henning Winker
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' mv <- ss3mvln(ss3sma, mc = 1000)
#' mv$refpts.iter
#' }
#' @export
ss3mvln <- function(ss3rep, mc = 1000, method = c("joint", "annual"), Fref = NULL, years = NULL, fwd = FALSE,
                    catch = c("Exp", "Obs", "kill_bio", "dead_bio", "ret_bio"), seed = 123, virgin = FALSE){
  method <- match.arg(method); catch <- match.arg(catch)
  if(catch == "dead_bio") catch <- "kill_bio"
  ss3rep$endyr <- .endyr(ss3rep)
  if(mc > 0 && is.null(ss3rep$CoVar))
    stop("CoVar required for MVLN draws: run SS3 with the hessian and read it with r4ss::SS_output(covar = TRUE)",
         call. = FALSE)
  if(method == "annual"){
    out <- .mvln_annual(ss3rep, Fref = Fref, years = years, mc = mc, forecast = fwd, virgin = virgin,
                        catch = catch, seed = seed)
    out$method <- "annual"
    return(out)
  }
  j <- .mvln_joint(ss3rep, Fref = Fref, years = years, mc = mc, forecast = fwd, catch = catch, seed = seed)
  b <- j$b
  rp <- function(d) FLCore::FLPar(Ftgt = d$Ftgt, Btgt = d$Btgt, MSY = d$MSY, B0 = d$B0, R0 = d$R0)
  rpi <- if(!is.null(j$rp)) FLCore::FLPar(array(t(as.matrix(j$rp[, .(Ftgt, Btgt, MSY, B0, R0)])), dim = c(5, mc),
                                                dimnames = list(params = c("Ftgt", "Btgt", "MSY", "B0", "R0"),
                                                                iter = seq_len(mc)))) else NULL
  structure(list(kb = j$kb, mle = j$mle, refpts = rp(j$rp.mle), refpts.iter = rpi,
                 basis = list(B = b$Bref, F = b$Fref, C = b$Cref, bbasis = b$bbasis, fbasis = b$fbasis, Btgtref = b$bref,
                              ssbunits = ss3rep$SpawnOutputUnits %||% "biomass"),
                 method = "joint"), class = "ss3mvln")
}

#' @export
print.ss3mvln <- function(x, ...){
  cat("SS3 MVLN (", x$method, "): ", if(is.null(x$kb)) 0 else max(x$kb$iter), " draws, years ", min(x$mle$year), "-",
      max(x$mle$year), "\nB ratio: ", x$basis$bbasis, " (B/", x$basis$B, "); F: ", x$basis$fbasis, " (F/", x$basis$F, ")\n",
      sep = "")
  print(x$refpts)
  invisible(x)
}


# year-by-year MVLN of ss3diags::SSdeltaMVLN (Walter and Winker 2019), without dplyr
.mvln_annual <- function(ss3rep, Fref = NULL, years = NULL, mc = 1000, forecast = FALSE, virgin = FALSE,
                         catch = "Exp", seed = 123){
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
  C <- .ss3_catch_year(ss3rep, catch)
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
                 basis = list(B = b$Bref, F = b$Fref, C = b$Cref, bbasis = b$bbasis, fbasis = b$fbasis, Btgtref = b$bref,
                              ssbunits = ss3rep$SpawnOutputUnits %||% "biomass")),
            class = "ss3mvln")
}


#' SS3 stock as FLStockR
#'
#' Builds the simplified (biomass-dynamic, age "1") \code{FLStockR} of \code{FLRef::ss2FLStockR()} from an SS3 report
#' or an \code{\link{ss3mvln}()} object: \code{ssb()} = SSB, \code{fbar()} = F, \code{rec()} = recruitment,
#' \code{catch()} = catch, with reference points (Ftgt, Btgt, MSY, B0, R0, and Blim, Bpa when given) in
#' \code{@refpts}, so SS3 runs plot and compare with JABBA3, FLicc and FLRef outputs.
#'
#' \itemize{
#'   \item \code{output = "mle"}: maximum likelihood estimates from the derived quantities (no CoVar needed), with
#'     log-scale SEs of SSB, F, recruitment, B/Btgt and F/Ftgt in \code{attr(x, "se.log")}; for runs without
#'     hessian/sdreport (derived quantities reported as 0) the time series are used (\code{source = "timeseries"}).
#'   \item \code{output = "iters"}: MVLN draws as iterations (\code{\link{ss3mvln}()}; \code{method = "joint"}: coherent
#'     trajectories and reference points with iterations; \code{"annual"}: SSdeltaMVLN, with the jointly sampled
#'     B/Btgt, F/Ftgt in \code{attr(x, "status")}); MLE in \code{attr(x, "mle")}.
#' }
#' With \code{rel = TRUE} the stock holds B/Btgt in place of SSB and F/Ftgt in place of F, with Btgt = Ftgt = 1 (as
#' \code{FLRef::ss2FLStockR(rel = TRUE)}).
#'
#' @param x output of \code{r4ss::SS_output()} or an \code{ss3mvln} object
#' @param output \code{"mle"} (default for an SS3 report) or \code{"iters"} (default for an ss3mvln object)
#' @param iters number of MVLN draws for \code{output = "iters"}
#' @param method MVLN method, \code{"joint"} (default) or \code{"annual"}
#' @param rel relative quantities (B/Btgt, F/Ftgt) in the stock
#' @param thin keep every \code{thin}-th iteration
#' @param Fref,fwd,catch passed to \code{\link{ss3mvln}()} (\code{fwd = TRUE}: add the SS3 forecast years)
#' @param blim,bpa optional Blim and Bpa as fractions of Btgt, added to the reference points
#' @param name stock name
#' @param source \code{"auto"}, \code{"derived"} or \code{"timeseries"} (MLE only)
#' @param btgt biomass target as a fraction of B0 for \code{source = "timeseries"} (default the SS3 \code{btarg})
#' @param area for area-structured models: return \code{FLStocks} with the stock (as \code{output = "mle"}) and one
#'   FLStockR per area (MLE). By area, SSB is the spawning output at spawning time, recruitment and catch (expected
#'   catch of the fleets of the area) are summed over seasons, and F is the F experienced in the area: the Fbar over the
#'   F ages of the area's fleets (\code{ss3_partialF(groups = "area", method = "fbar")}; apical F if the F basis is not
#'   age-based). This F does not add up to the stock F. Area reference points are approximations: B0 and R0 are the
#'   unfished SSB and recruitment of the area, Btgt (and Blim, Bpa) are scaled by the area's share of B0, Ftgt is
#'   the stock Ftgt, and MSY is not defined.
#' @return FLStockR
#' @author Henning Winker, Laurence T. Kell
#' @examples
#' data(ss3simple)
#' stk <- ss3stk(ss3simple)                    # MLE
#' sa <- ss3stk(readRDS(system.file("extdata", "pra.27.3a4a", "pra.27.3a4a.rds", package = "ss3diags2")),
#'              area = TRUE)                   # FLStocks: stock and areas
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' stk  <- ss3stk(ss3sma)                      # MLE with SEs
#' stkf <- ss3stk(ss3sma, output = "iters", fwd = TRUE)   # with forecast years
#' stki <- ss3stk(ss3sma, output = "iters")    # joint MVLN iterations
#' stki <- ss3stk(ss3mvln(ss3sma, method = "annual"))
#' }
#' @export
ss3stk <- function(x, output = NULL, iters = 1000, method = c("joint", "annual"), rel = FALSE, thin = 1, Fref = NULL,
                   fwd = FALSE, catch = c("Exp", "Obs", "kill_bio", "dead_bio", "ret_bio"), blim = NULL,
                   bpa = NULL, name = "ss3", source = c("auto", "derived", "timeseries"), btgt = NULL, area = FALSE){
  method <- match.arg(method); catch <- match.arg(catch); source <- match.arg(source)
  if(area){
    if(!.is_ss3rep(x)) stop("area = TRUE needs the SS3 report (r4ss::SS_output())", call. = FALSE)
    tot <- ss3stk(x, output = "mle", Fref = Fref, catch = catch, blim = blim, bpa = bpa, name = "stock",
                  source = source, btgt = btgt)
    return(.ss3stk_area(x, tot, blim, bpa))
  }
  if(catch == "dead_bio") catch <- "kill_bio"
  ss3rep <- NULL; se <- NULL
  if(inherits(x, "ss3mvln")){
    mv <- x; output <- output %||% (if(is.null(mv$kb)) "mle" else "iters")
  } else {
    ss3rep <- x; ss3rep$endyr <- .endyr(ss3rep); output <- output %||% "mle"
    output <- match.arg(output, c("mle", "iters"))
    if(output == "mle"){
      if(source == "auto") source <- if(.derived_ok(ss3rep)) "derived" else "timeseries"
      if(source == "timeseries"){
        message("derived quantities not available: trajectories from the time series; refpts B0, R0",
                if(isTRUE((btgt %||% ss3rep$btarg %||% -1) > 0)) " and Btgt" else "")
        mv <- .ss3_ts_mle(ss3rep, btgt); mv$method <- "timeseries"
      } else {
        j <- .mvln_joint(ss3rep, Fref = Fref, mc = 0, forecast = fwd, catch = catch)
        mv <- list(mle = j$mle, refpts = FLCore::FLPar(Ftgt = j$rp.mle$Ftgt, Btgt = j$rp.mle$Btgt, MSY = j$rp.mle$MSY,
                                                       B0 = j$rp.mle$B0, R0 = j$rp.mle$R0),
                   basis = list(B = j$b$Bref, F = j$b$Fref, C = j$b$Cref, bbasis = j$b$bbasis, fbasis = j$b$fbasis,
                                Btgtref = j$b$bref, ssbunits = ss3rep$SpawnOutputUnits %||% "biomass"), method = "mle")
        se <- .se_log(ss3rep, mv$mle$year)
      }
    } else mv <- ss3mvln(ss3rep, mc = iters, method = method, Fref = Fref, fwd = fwd, catch = catch)
  }
  output <- match.arg(output, c("mle", "iters"))
  if(output == "iters" && is.null(mv$kb)) stop("no draws in this ss3mvln object", call. = FALSE)
  d <- if(output == "iters") copy(mv$kb) else copy(mv$mle)[, iter := 1L]
  if(output == "iters" && thin > 1){
    keep <- unique(d$iter)[seq(1, length(unique(d$iter)), by = thin)]
    d <- d[iter %in% keep][, iter := match(iter, keep)]
  }
  if(rel) d[, `:=`(SSB = stock, F = harvest)]
  mk <- function(v, u) long2flq(d[, .(age = "1", year, iter, data = get(v))], quant = "age", units = u)
  stk <- .stk_build(d, name)
  # reference points: with iterations for the joint MVLN, else at the MLE
  rp <- if(output == "iters" && !is.null(mv$refpts.iter)){
    r <- mv$refpts.iter
    if(thin > 1) r <- FLCore::iter(r, keep)
    r
  } else mv$refpts
  if(rel){ rp["Btgt"] <- 1; rp["Ftgt"] <- 1 }
  extra <- c(Blim = blim %||% NA, Bpa = bpa %||% NA); extra <- extra[is.finite(extra)]
  if(length(extra)){
    nm <- c(dimnames(rp)$params, names(extra)); ni <- dim(rp)[length(dim(rp))]
    a <- rbind(matrix(rp@.Data, nrow = dim(rp)[1]), outer(extra, as.numeric(rp["Btgt"]))[, seq_len(ni), drop = FALSE])
    rp <- FLCore::FLPar(array(a, dim = c(length(nm), ni), dimnames = list(params = nm, iter = seq_len(ni))))
  }
  out <- FLCore::FLStockR(stk, refpts = rp)
  attr(out, "basis") <- mv$basis
  attr(out, "endyr") <- if(!is.null(ss3rep)) ss3rep$endyr else if("type" %in% names(mv$mle))
    max(mv$mle$year[mv$mle$type == "fit"]) else NULL
  if(output == "iters"){
    attr(out, "mle") <- mv$mle
    if(identical(mv$method, "annual"))       # SSdeltaMVLN: SSB and F draws marginal, keep the joint B/F ratios
      attr(out, "status") <- FLCore::FLQuants(stock = mk("stock", ""), harvest = mk("harvest", ""))
  }
  if(!is.null(se)) attr(out, "se.log") <- se
  if(!is.null(ss3rep)) attr(out, "eqcurve") <- tryCatch(ss3_eqcurve(ss3rep), error = function(e) NULL)
  out
}

# log-scale SEs of SSB, F, recruitment, B ratio from the derived quantities
.se_log <- function(ss3rep, yrs){
  hat <- as.data.table(ss3rep$derived_quants)
  lab <- .col(hat, c("Label", "LABEL")); val <- .col(hat, c("Value", "VALUE")); sdv <- .col(hat, c("StdDev", "STDDEV"))
  sel <- function(k){ i <- match(paste0(k, "_", yrs), lab); sqrt(log(1 + (sdv[i]/val[i])^2)) }
  sq <- function(v) long2flq(data.table(age = "1", year = yrs, data = v), quant = "age", units = "")
  FLCore::FLQuants(SSB = sq(sel("SSB")), F = sq(sel("F")), Rec = sq(sel("Recr")), stock = sq(sel("Bratio")),
                   harvest = sq(sel("F")))
}

#' Key quantities as FLQuants
#'
#' SSB, F, recruitment, catch and stock status (B/Btgt, F/Ftgt, B/B0, Catch/MSY) as \code{FLQuants} (with
#' iterations where the stock has them), for plotting with \code{ggplotFL} (\code{plot(ss3_flqs(x))}) or
#' \code{\link{ss3plot_trj}()}, as \code{JABBA3::jb3_flqs()}. Ratios are computed within each iteration, with
#' iteration-specific reference points for the joint MVLN.
#'
#' @param x FLStockR (e.g. \code{\link{ss3stk}()}), an SS3 report or an \code{ss3mvln} object (passed to \code{ss3stk()})
#' @param quants any of "SSB", "F", "Rec", "Catch", "B/Btgt", "F/Ftgt", "B/B0", "Catch/MSY"
#' @param ... passed to \code{ss3stk()} when \code{x} is not an FLStockR
#' @return FLQuants
#' @examples
#' data(ss3simple)
#' ss3_flqs(ss3stk(ss3simple), c("SSB", "F", "B/Btgt"))
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' ggplotFL::plot(ss3_flqs(ss3sma, output = "iters"))
#' }
#' @export
ss3_flqs <- function(x, quants = c("SSB", "F", "Rec", "Catch", "B/Btgt", "F/Ftgt"), ...){
  stk <- if(is(x, "FLStock")) x else ss3stk(x, ...)
  B <- FLCore::ssb(stk); F <- FLCore::fbar(stk); R <- FLCore::rec(stk); C <- FLCore::catch(stk)
  rp <- if(is(stk, "FLStockR")) FLCore::refpts(stk) else NULL
  rq <- function(q, p){
    if(is.null(rp) || !p %in% dimnames(rp)$params) return(q*NA_real_)
    r <- rp[p, ]
    if(dim(r)[length(dim(r))] > 1 && dim(q)[6] == dim(r)[length(dim(r))]) q %/% r else q/stats::median(as.numeric(r))
  }
  st <- attr(stk, "status")
  all <- list(SSB = B, F = F, Rec = R, Catch = C,
              `B/Btgt` = if(!is.null(st)) st$stock else rq(B, "Btgt"),
              `F/Ftgt` = if(!is.null(st)) st$harvest else rq(F, "Ftgt"),
              `B/B0` = rq(B, "B0"), `Catch/MSY` = rq(C, "MSY"))
  bad <- setdiff(quants, names(all)); if(length(bad)) stop("unknown quants: ", paste(bad, collapse = ", "), call. = FALSE)
  FLCore::FLQuants(all[quants])
}

utils::globalVariables(c("Ftgt", "Btgt", "MSY", "B0", "R0"))

# simplified (age "1") FLStock from a data.table year, iter, SSB, Recr, Catch, F
.stk_build <- function(d, name){
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
  stk
}

# stock by area: SSB, recruitment and catch by area, F experienced in each area (Fbar of the area fleets)
.ss3stk_area <- function(ss3rep, tot, blim = NULL, bpa = NULL){
  na <- ss3rep$nareas %||% 1
  if(na < 2) stop("the model has one area", call. = FALSE)
  yrs <- as.numeric(dimnames(tot)$year)
  ba <- ss3_byarea(ss3rep, years = yrs)
  long <- function(q, v) { x <- as.data.table(methods::as(ba[[q]], "data.frame"))[, .(year = as.numeric(as.character(year)), area = as.character(area), data)]
                           setnames(x, "data", v); x }
  d <- Reduce(function(a, b) merge(a, b, by = c("year", "area")), list(long("SSB", "SSB"), long("Rec", "Recr")))
  fl <- ss3_fleets(ss3rep)
  ct <- ss3_catch(ss3rep, what = "exp", long = TRUE)
  ct <- merge(ct[, .(c = sum(data)), by = .(name, year)], fl[, .(name, area = as.character(area))], by = "name")
  d <- merge(d, ct[, .(Catch = sum(c)), by = .(year, area)], by = c("year", "area"), all.x = TRUE)
  meth <- if(is.null(.fbar_ages(ss3rep))) "apical" else "fbar"
  fa <- ss3_partialF(ss3rep, groups = "area", method = meth, years = yrs, long = TRUE)
  d <- merge(d, fa[, .(year, area = sub("area ", "", name), F = data)], by = c("year", "area"), all.x = TRUE)
  d[is.na(Catch), Catch := 0][is.na(F), F := 0][, iter := 1L]
  # area reference points: unfished SSB and recruitment by area (virgin), Btgt in proportion to B0, stock Ftgt
  ts <- as.data.table(ss3rep$timeseries); sp <- ss3rep$spawnseas %||% 1
  v <- ts[Era == "VIRG"]
  b0a <- v[Seas == sp, .(b0 = sum(SpawnBio)), by = Area]; r0a <- v[, .(r0 = sum(Recruit_0)), by = Area]
  rp0 <- tot@refpts; bt <- .rp(tot, "Btgt"); B0 <- .rp(tot, "B0"); R0 <- .rp(tot, "R0")
  out <- lapply(sort(unique(d$area)), function(a){
    sh <- b0a[Area == as.numeric(a)]$b0/sum(b0a$b0); shr <- r0a[Area == as.numeric(a)]$r0/sum(r0a$r0)
    pv <- c(Ftgt = .rp(tot, "Ftgt"), Btgt = bt*sh, MSY = NA, B0 = B0*sh, R0 = R0*shr)
    ex <- c(Blim = blim %||% NA, Bpa = bpa %||% NA); ex <- ex[is.finite(ex)]*bt*sh
    pv <- c(pv, ex)
    st <- FLCore::FLStockR(.stk_build(d[area == a], paste("area", a)),
                           refpts = FLCore::FLPar(array(pv, dim = c(length(pv), 1),
                                                        dimnames = list(params = names(pv), iter = 1))))
    attr(st, "basis") <- c(attr(tot, "basis"), list(area = TRUE, Fmethod = meth))
    attr(st, "endyr") <- attr(tot, "endyr")
    st
  })
  names(out) <- paste("area", sort(unique(d$area)))
  FLCore::FLStocks(c(list(stock = tot), out))
}
