# lhpar.R - life-history parameters as FLPar

#' Life-history parameters
#'
#' Growth, maturity, weight-length and natural mortality parameters of an SS3 model as an \code{FLPar} with the
#' FLicc / FLRef names: \code{linf}, \code{k}, \code{t0} (von Bertalanffy; \code{t0} is the age at length 0 implied by
#' the SS3 L1-L2 parameterisation), \code{L50} and \code{mat.slope} (length-logistic maturity), \code{a}, \code{b}
#' (weight-length, kg and cm), \code{M} (natural mortality of the oldest ages), \code{CVmin}, \code{CVmax} (CV of
#' length at age), \code{A1}, \code{A2}, \code{L1}, \code{L2}. With two sexes the FLPar has a \code{unit} dimension
#' ("F", "M"), as the FLQuants of the compositions.
#'
#' Growth other than von Bertalanffy (SS3 Richards, age-specific K, growth cessation) is approximated by a
#' von Bertalanffy curve fitted to the SS3 mean length at age (\code{endgrowth}), with a message. Maturity other than
#' length-logistic leaves \code{L50} and \code{mat.slope} as NA.
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param year year of the (possibly time-varying) parameters (default: last year with parameters)
#' @return FLPar (params x unit with two sexes)
#' @examples
#' data(ss3ank)
#' ss3_lhpar(ss3ank)
#' @export
ss3_lhpar <- function(ss3rep, year = NULL){
  gp <- as.data.table(ss3rep$Growth_Parameters)
  if(!nrow(gp)) stop("no Growth_Parameters in ss3rep", call. = FALSE)
  yr <- year %||% max(gp$Yr)
  gp <- gp[Yr == max(gp$Yr[gp$Yr <= yr])]
  gp <- gp[!duplicated(Sex)]
  nsx <- ss3rep$nsexes %||% 1
  vb <- (ss3rep$GrowthModel_option %||% 1) == 1
  pars <- c("linf", "k", "t0", "L50", "mat.slope", "a", "b", "M", "CVmin", "CVmax", "A1", "A2", "L1", "L2")
  if(nsx == 2 && !2 %in% gp$Sex){
    message("no male growth parameters in Growth_Parameters: male values set to the female values")
    gp <- rbind(gp, copy(gp[Sex == 1])[, Sex := 2])
  }
  one <- function(s){
    g <- gp[Sex == s]
    linf <- g$Linf; k <- g$K; t0 <- .col(g, "A_a_L0", NA_real_)
    if(!vb){
      eg <- as.data.table(ss3rep$endgrowth)[Sex == s & Seas == 1]
      eg <- eg[!duplicated(Age_Beg)]
      f <- stats::nls(Len_Mid ~ li*(1 - exp(-kk*(Age_Mid - tt))), data = eg,
                      start = list(li = max(eg$Len_Mid), kk = 0.2, tt = 0))
      cf <- stats::coef(f); linf <- cf[["li"]]; k <- cf[["kk"]]; t0 <- cf[["tt"]]
    }
    fem <- gp[Sex == 1]
    c(linf = linf, k = k, t0 = t0,
      L50 = if((ss3rep$Maturity_option %||% 1) == 1) fem$Mat1 else NA_real_,
      mat.slope = if((ss3rep$Maturity_option %||% 1) == 1) fem$Mat2 else NA_real_,
      a = g$WtLen1, b = g$WtLen2, M = .col(g, "M_nages", NA_real_), CVmin = g$CVmin, CVmax = g$CVmax,
      A1 = g$A1, A2 = g$A2, L1 = g$L_a_A1, L2 = g$L_a_A2)
  }
  if(!vb) message("growth is not von Bertalanffy: linf, k, t0 fitted to the SS3 mean length at age")
  if(nsx == 1) return(FLCore::FLPar(one(1)[pars], params = pars))
  arr <- array(c(one(1)[pars], one(2)[pars]), dim = c(length(pars), 2, 1),
               dimnames = list(params = pars, unit = c("F", "M"), iter = 1))
  FLCore::FLPar(arr)
}

# von Bertalanffy parameters for unit u from an FLPar with or without a unit dimension
.vb_unit <- function(lh, u){
  u <- as.character(u)
  d <- dimnames(lh)
  get <- function(p){
    a <- lh@.Data
    if("unit" %in% names(d)){
      uu <- if(isTRUE(u %in% d$unit)) u else d$unit[1]
      as.numeric(a[p, uu, 1])
    } else as.numeric(a[p, 1])
  }
  c(linf = get("linf"), k = get("k"), t0 = if("t0" %in% d$params) get("t0") else 0)
}

utils::globalVariables(c("Yr", "Sex", "Seas", "Age_Beg", "Age_Mid", "Len_Mid"))
