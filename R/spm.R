# spm.R - Pella-Tomlinson surplus production function fitted to the SS3 equilibrium yield curve

# Pella-Tomlinson surplus production, with the Fox limit at m = 1
.pt <- function(B, r, K, m){
  sp <- if(abs(m - 1) < 1e-6) r*B*log(K/pmax(B, 1e-12)) else r/(m - 1)*B*(1 - (B/K)^(m - 1))
  sp[B <= 0] <- 0
  sp
}

# shape m from Bmsy/K: Bmsy/K = m^(-1/(m-1)) (1/e at m = 1, Fox)
.m_from_bk <- function(bk){
  if(abs(bk - exp(-1)) < 1e-8) return(1)
  stats::uniroot(function(m) (if(abs(m - 1) < 1e-8) exp(-1) else m^(-1/(m - 1))) - bk, c(1e-4, 100), tol = 1e-10)$root
}

.pt_refs <- function(r, K, m){
  bk <- if(abs(m - 1) < 1e-6) exp(-1) else m^(-1/(m - 1))
  msy <- .pt(bk*K, r, K, m)
  c(MSY = msy, Bmsy = bk*K, BmsyK = bk, Hmsy = msy/(bk*K))
}

#' Surplus production function from the SS3 equilibrium curve
#'
#' Fits a Pella-Tomlinson production function, \eqn{SP(B) = r/(m-1) B (1 - (B/K)^{m-1})} (Fox at m = 1), to the SS3
#' equilibrium yield curve (yield against SSB, \code{\link{ss3_eqcurve}()}), as \code{FLRef::asem2spm()} does for an
#' FLBRP. The fit summarises the age-structured equilibrium in three parameters and extrapolates it beyond the range of
#' F profiled by SS3 (down to SSB = 0 and up to K). It also gives SS3-based priors for r, K and the shape m of
#' JABBA3 (with \code{model = "SB"}-type runs, as the curve is in spawning biomass).
#'
#' Extrapolation below the profiled range rests on the Pella-Tomlinson shape: with m < 1 (Bmsy/K < 0.37) production per
#' unit biomass rises without bound as SSB goes to 0, so the curve is optimistic at very low SSB, where the SS3 curve
#' falls to zero at the crash F.
#'
#' With \code{K = "free"} (default) K is estimated with r and m; with \code{K = "B0"} it is fixed at the SS3 unfished
#' SSB. The fit is checked against the SS3 benchmarks (MSY, SSB_MSY, B0): a relative difference above \code{tol}
#' gives a warning, as it means the equilibrium table and the SS3 MSY calculation use different settings (e.g.
#' benchmark years for biology or selectivity) or that the curve is not of Pella-Tomlinson shape.
#'
#' @param x output of \code{r4ss::SS_output()} or of \code{\link{ss3_eqcurve}()}
#' @param K \code{"free"} (default) or \code{"B0"}
#' @param tol tolerance for the check against the SS3 benchmarks
#' @param cv optional CVs for the priors; default from the SS3 standard errors of B0 (K) and of MSY and B0 (r)
#' @return list of class \code{ss3spm}: \code{pars} (FLPar: r, K, m, MSY, Bmsy, BmsyK, Hmsy = MSY/Bmsy),
#'   \code{curve} (data.table SSB, yield on 0..K), \code{eq} (the SS3 curve with fitted values), \code{check}
#'   (SS3 vs SPM benchmarks), \code{r2}, and \code{priors} (list of \code{c(median, log.sd, on)} for r, K and m, the
#'   format of \code{JABBA3::jb3_priors()}; m is off by default; r and K are off when SSB is in numbers)
#' @examples
#' \dontrun{
#' data(pac.hke, package = "ss3diags")
#' sp <- ss3_spm(ss3phk)
#' sp
#' # JABBA3::jb3_priors(data, r = sp$priors$r, K = sp$priors$K, m = sp$priors$m)
#' }
#' @export
ss3_spm <- function(x, K = c("free", "B0"), tol = 0.05, cv = NULL){
  K <- match.arg(K)
  eq <- if(.is_ss3rep(x)) ss3_eqcurve(x) else x
  if(!all(c("SSB", "yield") %in% names(eq))) stop("x must be an ss3rep or the output of ss3_eqcurve()", call. = FALSE)
  d <- as.data.table(eq)[SSB >= 0 & yield >= 0]
  ss3 <- c(MSY = attr(eq, "MSY") %||% NA, Bmsy = attr(eq, "Bmsy") %||% NA, B0 = attr(eq, "B0") %||% NA)
  # starting values: K at the unfished end of the curve, shape from the curve maximum
  K0 <- max(d$SSB); i <- which.max(d$yield)
  m0 <- .m_from_bk(min(max(d$SSB[i]/K0, 0.05), 0.95))
  r0 <- d$yield[i]/(K0*.pt_refs(1, 1, m0)[["MSY"]])
  sc <- max(d$yield)
  obj <- function(p){
    kk <- if(K == "free") exp(p[3]) else ss3[["B0"]]
    sum(((d$yield - .pt(d$SSB, exp(p[1]), kk, exp(p[2])))/sc)^2)
  }
  p0 <- if(K == "free") log(c(r0, m0, K0)) else log(c(r0, m0))
  fit <- stats::optim(p0, obj, method = "Nelder-Mead", control = list(maxit = 5000, reltol = 1e-12))
  fit <- stats::optim(fit$par, obj, method = "BFGS", control = list(maxit = 1000, reltol = 1e-14))
  r <- exp(fit$par[1]); m <- exp(fit$par[2]); kk <- if(K == "free") exp(fit$par[3]) else ss3[["B0"]]
  rf <- .pt_refs(r, kk, m)
  d[, fit := .pt(SSB, r, kk, m)]
  r2 <- 1 - sum((d$yield - d$fit)^2)/sum((d$yield - mean(d$yield))^2)
  check <- data.table(quantity = c("MSY", "Bmsy", "B0", "Bmsy/B0"),
                      SS3 = c(ss3[["MSY"]], ss3[["Bmsy"]], ss3[["B0"]], ss3[["Bmsy"]]/ss3[["B0"]]),
                      SPM = c(rf[["MSY"]], rf[["Bmsy"]], kk, rf[["BmsyK"]]))
  check[, rel.diff := SPM/SS3 - 1]
  bad <- check[is.finite(rel.diff) & abs(rel.diff) > tol]
  if(nrow(bad)) warning("SPM fit and SS3 benchmarks differ by more than ", 100*tol, "% for ",
                        paste(bad$quantity, collapse = ", "), ": check the SS3 equilibrium settings (see $check)",
                        call. = FALSE)
  # priors (lognormal): CVs from SS3 standard errors unless given
  sdq <- attr(eq, "sd") %||% c(MSY = NA, B0 = NA)
  cvK <- cv[["K"]] %||% (sdq[["B0"]]/ss3[["B0"]]); cvR <- cv[["r"]] %||% sqrt((sdq[["MSY"]]/ss3[["MSY"]])^2 + cvK^2)
  lsd <- function(cv) if(is.finite(cv)) sqrt(log(1 + cv^2)) else 0.3
  priors <- list(r = c(median = r, log.sd = lsd(cvR), on = 1), K = c(median = kk, log.sd = lsd(cvK), on = 1),
                 m = c(median = m, log.sd = lsd(cv[["m"]] %||% 0.1), on = 0))
  if(identical(attr(eq, "ssbunits"), "numbers")){
    warning("SSB is spawning output in numbers: r and K mix units (yield in weight), so only the shape (m, BmsyK) ",
            "carries over; r and K priors are switched off", call. = FALSE)
    priors$r[["on"]] <- 0; priors$K[["on"]] <- 0
  }
  bs <- seq(0, kk, length.out = 400)
  structure(list(pars = FLCore::FLPar(r = r, K = kk, m = m, MSY = rf[["MSY"]], Bmsy = rf[["Bmsy"]],
                                      BmsyK = rf[["BmsyK"]], Hmsy = rf[["Hmsy"]]),
                 curve = data.table(SSB = bs, yield = .pt(bs, r, kk, m)), eq = d[], check = check[], r2 = r2,
                 priors = priors, K = K, ssbunits = attr(eq, "ssbunits")),
            class = "ss3spm")
}

#' @export
print.ss3spm <- function(x, ...){
  p <- as.numeric(x$pars); names(p) <- dimnames(x$pars)$params
  cat("Pella-Tomlinson fit to the SS3 equilibrium curve (K ", x$K, "), R2 = ", format(x$r2, digits = 5), "\n", sep = "")
  print(signif(p, 4))
  cat("\nCheck against SS3 benchmarks:\n"); print(x$check, digits = 4)
  invisible(x)
}

utils::globalVariables(c("rel.diff", "loop", "f0", "f1"))
