# lhpar.R - life-history parameters as FLPar

#' Life-history parameters
#'
#' Growth, maturity, weight-length and natural mortality parameters of an SS3 model as an \code{FLPar} with the
#' FLicc / FLRef names: \code{linf}, \code{k}, \code{t0} (von Bertalanffy; \code{t0} is the age at length 0 implied by
#' the SS3 L1-L2 parameterisation), \code{L50} and \code{mat.slope} (length-logistic maturity), \code{a}, \code{b}
#' (weight-length, kg and cm), \code{M} (natural mortality of the oldest ages), \code{CVmin}, \code{CVmax} (CV of
#' length at age), \code{A1}, \code{A2}, \code{L1}, \code{L2}, \code{K} and \code{richards} (SS3 growth in the Schnute
#' (1981) parameterisation; \code{richards} = 1 for von Bertalanffy, the Richards shape for SS3 growth option 2, NA
#' otherwise), \code{settle}, the time of settlement as a fraction of the year ((Settle_Month - 1)/12),
#' \code{cvtype} (SS3 CV pattern: 0 CV = f(length at age), 1 CV = f(age), 2 SD = f(length at age), 3 SD = f(age);
#' \code{CVmin}, \code{CVmax} are SDs for 2 and 3) and \code{amax} (accumulator age, used for A2 = 999). With two
#' sexes the FLPar has a \code{unit} dimension ("F", "M"), as the FLQuants of the compositions.
#'
#' Ages: SS3 data and \code{endgrowth} use calendar ages (age 0 from 1 January), and so do \code{linf}, \code{k},
#' \code{t0} (SS3 \code{A_a_L0}); the Schnute parameters \code{A1} and \code{A2} are ages since settlement, so the
#' growth curves evaluate them at calendar age - \code{settle}. This matters for seasonal models with settlement after
#' January (e.g. Settle_Month 6: fish are 5 months younger than their calendar age).
#'
#' For growth other than von Bertalanffy, \code{linf}, \code{k}, \code{t0} are a von Bertalanffy curve fitted to the
#' SS3 mean length at age (\code{endgrowth}), with a message; the growth curves of the plots use the SS3 Richards
#' (Schnute) curve when \code{richards} is available. Maturity other than
#' length-logistic leaves \code{L50} and \code{mat.slope} as NA.
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param year year of the (possibly time-varying) parameters (default: last year with parameters)
#' @return FLPar (params x unit with two sexes)
#' @examples
#' data(ss3simple)
#' ss3_lhpar(ss3simple)
#' hke <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18.rds", package = "ss3diags2"))
#' ss3_lhpar(hke)
#' @export
ss3_lhpar <- function(ss3rep, year = NULL){
  gp <- as.data.table(ss3rep$Growth_Parameters)
  if(!nrow(gp)) stop("no Growth_Parameters in ss3rep", call. = FALSE)
  yr <- year %||% max(gp$Yr)
  gp <- gp[Yr == max(gp$Yr[gp$Yr <= yr])]
  gp <- gp[!duplicated(Sex)]
  nsx <- ss3rep$nsexes %||% 1
  vb <- (ss3rep$GrowthModel_option %||% 1) == 1
  gm <- ss3rep$GrowthModel_option %||% 1
  pars <- c("linf", "k", "t0", "L50", "mat.slope", "a", "b", "M", "CVmin", "CVmax", "A1", "A2", "L1", "L2", "K",
            "richards", "settle", "cvtype", "amax")
  # CV of length at age: 0 CV = f(LAA), 1 CV = f(A), 2 SD = f(LAA), 3 SD = f(A) (SS3 CV_Growth_Pattern); amax is the
  # accumulator age, which SS3 uses for A2 = 999 (A2 = Linf)
  cvt <- tolower(ss3rep$growthCVtype %||% "cv=f(laa)")
  cvtype <- if(grepl("^sd", cvt)) 2 else 0
  if(grepl("\\(a\\)$", cvt)) cvtype <- cvtype + 1
  amax <- as.numeric(ss3rep$accuage %||% NA)
  # settlement: SS3 ages are calendar ages (age 0 from 1 January) but growth starts at settlement, so A1 and A2 are
  # ages since settlement; settle = (Settle_Month - 1)/12 of the main settlement event (0 if not reported)
  rd <- ss3rep$recruitment_dist$recruit_dist
  settle <- 0
  if(is.data.frame(rd) && "Settle_Month" %in% names(rd) && nrow(rd)){
    w <- rowSums(rd[, grep("^recr_dist", names(rd)), drop = FALSE], na.rm = TRUE)   # NA for males of 1-sex models
    i <- if(any(w > 0)) which.max(w) else 1
    settle <- (suppressWarnings(as.numeric(as.character(rd$Settle_Month[i]))) - 1)/12
    if(length(settle) != 1 || !is.finite(settle) || settle < 0) settle <- 0
  }
  # SS3 growth (Schnute 1981 parameterisation): K and the Richards shape (1 = von Bertalanffy; NA for other models)
  rich <- function(s){
    if(gm == 1) return(1)
    if(gm != 2) return(NA_real_)
    pp <- ss3rep$parameters
    lab <- grep(paste0("^Richards_", c("Fem", "Mal")[s]), rownames(pp), value = TRUE)
    if(length(lab)) pp[lab[1], "Value"] else NA_real_
  }
  if(nsx == 2 && !2 %in% gp$Sex){
    # older SS3 reports list female growth only: male growth from the parameters, as direct values or as exp offsets
    # from the female values, whichever reproduces the male mean length at age of endgrowth
    gm2 <- copy(gp[Sex == 1])[, Sex := 2]
    pp <- ss3rep$parameters
    pv <- function(lab) { l <- grep(lab, rownames(pp), value = TRUE); if(length(l)) pp[l[1], "Value"] else NA_real_ }
    mal <- c(L_a_A1 = pv("^L_at_Amin_Mal"), L_a_A2 = pv("^L_at_Amax_Mal"), K = pv("^VonBert_K_Mal"),
             CVmin = pv("^CV_young_Mal"), CVmax = pv("^CV_old_Mal"))
    fem <- unlist(gp[Sex == 1, names(mal), with = FALSE])
    eg <- as.data.table(ss3rep$endgrowth)[Sex == 2 & Seas == 1]
    vbL <- function(v) v[["L_a_A2"]] + (v[["L_a_A1"]] - v[["L_a_A2"]])*exp(-v[["K"]]*(eg$Age_Mid - gm2$A1))
    if(all(is.finite(mal[1:3])) && nrow(eg) && gm == 1){
      cand <- list(direct = mal, offset = fem*exp(mal))
      err <- vapply(cand, function(v) mean(abs(vbL(v) - eg$Len_Mid)), 0)
      v <- cand[[which.min(err)]]
      for(n in names(v)) if(is.finite(v[[n]])) gm2[, (n) := v[[n]]]
      gm2[, Linf := L_a_A2]
      gm2[, A_a_L0 := A1 + settle + log(1 - L_a_A1/L_a_A2)/K]          # calendar age, as SS3
      message("male growth from the parameters (", names(cand)[which.min(err)], " values): not listed in Growth_Parameters")
    } else message("no male growth parameters in Growth_Parameters: male values set to the female values")
    gp <- rbind(gp, gm2, fill = TRUE)
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
      A1 = g$A1, A2 = g$A2, L1 = g$L_a_A1, L2 = g$L_a_A2, K = g$K, richards = rich(s), settle = settle,
      cvtype = cvtype, amax = amax)
  }
  if(!vb) message(if(gm == 2) "Richards growth: linf, k, t0 are a von Bertalanffy fit to the SS3 mean length at age; growth curves use the Richards (Schnute) parameters A1, A2, L1, L2, K, richards"
                  else "growth is not von Bertalanffy: linf, k, t0 fitted to the SS3 mean length at age")
  if(nsx == 1) return(FLCore::FLPar(one(1)[pars], params = pars))
  arr <- array(c(one(1)[pars], one(2)[pars]), dim = c(length(pars), 2, 1),
               dimnames = list(params = pars, unit = c("F", "M"), iter = 1))
  FLCore::FLPar(arr)
}

#' Von Bertalanffy parameters of a unit
#'
#' Von Bertalanffy parameters for unit u from an FLPar with or without a unit dimension
#'
#' @param lh FLPar of life-history parameters (\code{ss3_lhpar()})
#' @param u unit (\code{"F"}, \code{"M"} or \code{"unique"})
#' @return named numeric linf, k, t0
#' @noRd
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
