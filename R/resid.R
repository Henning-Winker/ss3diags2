# resid.R - Pearson and one-step-ahead (OSA) residuals of the SS3 compositions

#' Beta-binomial cdf (as FLicc)
#'
#' Beta-binomial cdf at x for N trials, shapes a, b (as FLicc)
#'
#' @param x counts
#' @param N number of trials
#' @param a shape a
#' @param b shape b
#' @return probabilities
#' @noRd
.pbetabinom <- function(x, N, a, b){
  if(N <= 0) return(rep(1, length(x)))
  k <- 0:N
  lp <- lchoose(N, k) + lbeta(k + a, N - k + b) - lbeta(a, b)
  cdf <- cumsum(exp(lp - max(lp))); cdf <- cdf/cdf[length(cdf)]
  ifelse(x < 0, 0, ifelse(x >= N, 1, cdf[pmin(floor(x), N) + 1]))
}

#' OSA residuals of one composition
#'
#' OSA residuals of one composition: counts o, expected proportions p, DM precision a0 (Inf = multinomial)
#'
#' @param o observed counts
#' @param p expected proportions
#' @param a0 Dirichlet-multinomial precision (Inf = multinomial)
#' @return numeric residuals
#' @noRd
.osa_comp <- function(o, p, a0 = Inf){
  K <- length(o); r <- rep(NA_real_, K)
  p <- pmax(p, 1e-12); p <- p/sum(p)
  N <- sum(o); rem_n <- N; rem_p <- 1
  for(k in seq_len(K - 1)){
    if(rem_n <= 0 || rem_p <= 1e-12) break
    if(is.finite(a0)){
      a <- a0*p[k]; b <- a0*max(rem_p - p[k], 1e-12)
      lo <- .pbetabinom(o[k] - 1, rem_n, a, b); hi <- .pbetabinom(o[k], rem_n, a, b)
    } else {
      q <- min(max(p[k]/rem_p, 1e-12), 1 - 1e-12)
      lo <- stats::pbinom(o[k] - 1, rem_n, q); hi <- stats::pbinom(o[k], rem_n, q)
    }
    u <- lo + stats::runif(1)*(hi - lo)
    r[k] <- stats::qnorm(min(max(u, 1e-12), 1 - 1e-12))
    rem_n <- rem_n - o[k]; rem_p <- rem_p - p[k]
  }
  r
}

#' Residuals of the SS3 length or age compositions
#'
#' Pearson or one-step-ahead (OSA) residuals of the compositions, consistent with the SS3 likelihood: multinomial with
#' the adjusted sample size, or Dirichlet-multinomial with the precision implied by the SS3 effective sample size
#' (\eqn{\alpha_0 = N (N_{eff} - 1)/(N - N_{eff})} from \code{Nsamp_DM} / \code{DM_effN}, valid for both SS3 DM
#' parameterisations).
#'
#' OSA residuals (Trijoulet et al. 2023, Fish. Res. 257: 106487) decompose each composition along increasing bins
#' (females before males in joint-sex compositions) into conditional binomials (multinomial) or beta-binomials
#' (Dirichlet-multinomial) and return randomised quantile residuals \eqn{\Phi^{-1}(u)}, independent standard normal
#' under a correct model; the last bin of each composition has no residual. Counts are \code{round(obs x N)}. Pearson
#' residuals are those reported by SS3 (\code{Pearson}), correlated within a composition and best read as bubbles.
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param quant \code{"len"} (default), \code{"age"}, \code{"con"} (conditional age-at-length: each length bin of a
#'   sample is its own multinomial over ages, so OSA residuals are computed per length bin; in the FLQuants the bins of
#'   an age x year cell are combined (Stouffer), \code{long = TRUE} keeps the length bin, for
#'   \code{\link{ss3plot_caal}()}) or \code{"size"} (generalised size compositions)
#' @param type \code{"osa"} (default) or \code{"pearson"}
#' @param fleets fleet numbers or names (default all)
#' @param min.exp OSA residuals of bins with zero observations and an expected count (proportion x N) below
#'   \code{min.exp} are set to NA: the randomised residual of a bin that cannot be observed is pure noise (e.g. old
#'   ages in small length bins of conditional age-at-length data, empty tails of length compositions); 0 keeps all
#' @param nsamp sample size of the compositions: \code{"adj"} (default, the adjusted sample size of the SS3
#'   likelihood) or \code{"input"} (the input sample size, e.g. the number of aged fish per length bin). With strong
#'   down-weighting (variance adjustments well below 1) the adjusted sample size can round to zero fish, leaving few
#'   OSA residuals; \code{"input"} shows the residuals of the data as sampled
#' @param seed random seed for the OSA randomisation
#' @param long return the long data.table instead of FLQuants
#' @return FLQuants of residuals (one per fleet and partition, as \code{\link{ss3_lencomp}()}: len/age x year x unit x
#'   season x area), or a data.table with the FLQuant dimension columns, \code{data} (residual), \code{obs}, \code{exp}
#'   (counts), \code{n}, \code{a0}, \code{fleet}, \code{name}, \code{partname}; attribute \code{type}
#' @examples
#' data(ss3simple)
#' r <- ss3_resid(ss3simple)
#' names(r)
#' ss3plot_resid(r)
#' ss3_resid(ss3simple, "age", type = "pearson")
#' @export
ss3_resid <- function(ss3rep, quant = c("len", "age", "con", "size"), type = c("osa", "pearson"), fleets = NULL,
                      seed = 1, min.exp = 0.01, nsamp = c("adj", "input"), long = FALSE){
  nsamp <- match.arg(nsamp)
  quant <- match.arg(quant); type <- match.arg(type)
  db <- .comp_db(ss3rep, quant)
  qn <- if(quant == "con") "age" else quant
  fleets <- .pick_fleets(ss3rep, fleets, .fleetnames(ss3rep)[sort(unique(db$Fleet))])
  db <- db[.fleetnames(ss3rep)[Fleet] %in% fleets]
  ns <- ss3rep$nseasons %||% 1; na <- ss3rep$nareas %||% 1; nsx <- ss3rep$nsexes %||% 1
  dm <- .col(db, c("Nsamp_DM", "DM_effN"))
  d <- data.table(fleet = db$Fleet, name = .fleetnames(ss3rep)[db$Fleet], q = db$Bin, year = db$Yr,
                  unit = .unit_lab(.col(db, "Sex", 1), .col(db, "Sexes", 0), nsx), sex = .col(db, "Sex", 1),
                  season = .season_lab(db$Seas, ns), area = .area_lab(.col(db, "Area", 1), na), iter = "1",
                  part = .col(db, "Part", 0), month = .col(db, "Month", 0), sexes = .col(db, "Sexes", 0),
                  ageerr = .col(db, "Ageerr", 0), repl = .col(db, "Repl.", 1), lbin = .col(db, "Lbin_lo", NA_real_),
                  method = .col(db, "method", 1),
                  pobs = db$Obs, pexp = db$Exp, pearson = .col(db, "Pearson", NA_real_), N = if(nsamp == "adj") db$Nsamp_adj else .col(db, "Nsamp_in", db$Nsamp_adj), Nadj = db$Nsamp_adj,
                  neff = dm %||% NA_real_)
  d[, partname := name]
  d[, np := uniqueN(part), by = name][np > 1, partname := paste0(name, ".", .part_lab[as.character(part)])][, np := NULL]
  # one composition = fleet x year x season x month x partition x sex option x ageing error x replicate
  # (conditional age-at-length: one composition per length bin; size: per size method)
  d[, comp := .GRP, by = .(fleet, year, season, area, month, part, sexes, ageerr, repl, lbin, method)]
  setorder(d, comp, sex, q)
  d[, a0 := { ne <- neff[1]; n <- Nadj[1]
              if(is.finite(ne) && ne > 1 && ne < n) n*(ne - 1)/(n - ne) else Inf }, by = comp]
  d[, `:=`(obs = round(pobs*N), exp = pexp*N, n = N)]
  if(type == "osa"){
    old <- if(exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if(!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
    d[, data := .osa_comp(obs, pexp, a0[1]), by = comp]
    # cells with nothing observed and (almost) nothing expected carry no information: their randomised residual
    # is pure noise (u ~ U(0, 1)), so it is dropped
    if(min.exp > 0) d[obs == 0 & exp < min.exp, data := NA_real_]
  } else d[, data := pearson]
  setnames(d, "q", qn)
  if(quant == "con" && long){
    d <- d[, c("fleet", "name", "partname", "lbin", qn, .flqdims, "data", "obs", "exp", "n", "a0", "pobs", "pexp"),
           with = FALSE]
    setattr(d, "type", type); setattr(d, "quant", "con")
    return(d[])
  }
  quant <- qn
  # several compositions in one FLQuant cell (months within a season, replicates, length bins of conditional
  # age-at-length): Stouffer combination
  by <- c("partname", quant, .flqdims)
  if(anyDuplicated(d[is.finite(data)], by = by))
    d <- d[, .(fleet = fleet[1], name = name[1], data = if(all(is.na(data))) NA_real_ else
                 sum(data, na.rm = TRUE)/sqrt(sum(is.finite(data))),
               obs = sum(obs), exp = sum(exp), n = sum(n[!duplicated(comp)]), a0 = a0[1]), by = by]
  d <- d[, c("fleet", "name", "partname", quant, .flqdims, "data", "obs", "exp", "n", "a0"), with = FALSE]
  setattr(d, "type", type)
  if(long) return(d[])
  # compositions without any informative residual (e.g. adjusted sample sizes of 1): dropped with a hint
  ok <- unique(d[is.finite(data)]$partname); bad <- setdiff(unique(d$partname), ok)
  if(length(bad)){
    hint <- if(nsamp == "adj" && type == "osa" && all(d$n <= 1, na.rm = TRUE))
      " (adjusted sample sizes of 1: try nsamp = \"input\")" else ""
    if(!length(ok)) stop("no informative residuals", hint, call. = FALSE)
    warning("no informative residuals for ", paste(bad, collapse = ", "), hint, call. = FALSE)
  }
  out <- long2flqs(d[is.finite(data)], quant = quant, units = paste(type, "residual"), by = "partname",
                   order = intersect(unique(d$partname), ok))
  attr(out, "type") <- type
  out
}

utils::globalVariables(c("comp", "a0", "neff", "N", "pobs", "repl", "np"))

#' Plot residuals of conditional age-at-length data
#'
#' OSA (or Pearson) residuals of conditional age-at-length compositions by age (x) and length bin (y), one panel per
#' year, for one fleet: bubbles sized by the absolute residual, blue for observed > expected, red for observed <
#' expected. Growth curves (\code{lhpar}) show the mean length at age with \eqn{\pm} 2 SD bands; diamonds and bars
#' give the observed and expected mean length at age of the aged fish.
#'
#' @param x long data.table from \code{ss3_resid(ss3rep, quant = "con", long = TRUE)}
#' @param fleet fleet name (default: the first fleet with data)
#' @param years years to plot (default: all with data)
#' @param units units (sexes) to plot (default all; two sexes are drawn side by side at each age, females as circles,
#'   males as triangles)
#' @param lhpar optional FLPar from \code{\link{ss3_lhpar}()} (SS3 von Bertalanffy or Richards growth) or the SS3 report
#'   (SS3 mean length at age) for the growth curves by unit (solid females, dashed males)
#' @param timing fraction of the year at sampling: the curve is drawn at age class + \code{timing} (default mid-year)
#' @param band shade \eqn{\pm} 2 SD of length at age around the curves (SS3 CV of length at age)
#' @param means add the mean length at age of the aged fish per year: observed (diamonds) and expected (bars), both
#'   under the sampling by length bin of the data (so they compare with each other, not with the population curve)
#' @param cap residuals capped at this absolute value for the bubble size
#' @param ncol number of panel columns
#' @param text show the title, subtitle and caption (\code{FALSE}: none; labels inside the plot and legends stay)
#' @return ggplot
#' @examples
#' hke <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18.rds", package = "ss3diags2"))
#' cl <- ss3_resid(hke, quant = "con", long = TRUE, nsamp = "input")
#' ss3plot_caal(cl, fleet = "ITA_OTB_18", years = 2018:2023, lhpar = ss3_lhpar(hke))
#' @export
ss3plot_caal <- function(x, fleet = NULL, years = NULL, units = NULL, lhpar = NULL, timing = 0.5, band = TRUE,
                         means = TRUE, cap = 4, ncol = NULL, text = TRUE){
  .drop_text((function(){
  if(!is.data.frame(x) || !"lbin" %in% names(x))
    stop("x must be the long table of ss3_resid(ss3rep, quant = \"con\", long = TRUE)", call. = FALSE)
  d0 <- as.data.table(x)
  fl <- fleet %||% d0$partname[1]
  d0 <- d0[partname %in% fl | name %in% fl]
  if(!is.null(years)) d0 <- d0[year %in% years]
  if(!is.null(units)) d0 <- d0[unit %in% units]
  if(!nrow(d0)) stop("no residuals for this fleet / years", call. = FALSE)
  # length bins at their mid points
  lb <- sort(unique(d0$lbin)); w <- if(length(lb) > 1) stats::median(diff(lb)) else 0
  d0[, `:=`(age = as.numeric(as.character(age)), L = lbin + w/2)]
  us <- sort(unique(d0$unit)); nu <- length(us)
  off <- function(u) if(nu > 1) (match(u, us) - (nu + 1)/2)*0.3 else 0
  d <- d0[is.finite(data)]
  d[, `:=`(r = pmin(abs(data), cap), sign = ifelse(data > 0, "obs > exp", "obs < exp"), xa = age + off(unit))]
  shp <- stats::setNames(c(21, 24, 22)[seq_len(nu)], us)
  p <- ggplot2::ggplot(d, ggplot2::aes(xa, L))
  if(!is.null(lhpar)){
    # mean length of an age class at the sampling time: curve at real age = age class + timing; band +/- 2 SD of
    # length at age (SS3 CV of length at age)
    gc <- .growth_curves(lhpar, maxage = max(d0$age) + 1)[, age := age - timing][age >= -0.5]
    gc <- gc[L <= max(d0$L)*1.1]
    if(band && "sd" %in% names(gc) && any(is.finite(gc$sd)))
      p <- p + ggplot2::geom_ribbon(data = gc[is.finite(sd)], ggplot2::aes(age, ymin = pmax(L - 2*sd, 0), ymax = L + 2*sd,
                                                                         group = unit), fill = "grey50", alpha = 0.12,
                                    inherit.aes = FALSE)
    p <- p + ggplot2::geom_line(data = gc, ggplot2::aes(age, L, linetype = unit), colour = "grey40", inherit.aes = FALSE)
  }
  p <- p +
    ggplot2::geom_point(ggplot2::aes(size = r, fill = sign, shape = unit), colour = "grey30", stroke = 0.2, alpha = 0.8) +
    ggplot2::scale_shape_manual(values = shp, guide = if(nu > 1) ggplot2::guide_legend(override.aes = list(size = 3, fill = "grey70")) else "none") +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(shape = 21, size = 3))) +
    ggplot2::scale_fill_manual(values = c(`obs > exp` = "#3C78B4", `obs < exp` = "#D2604B"), name = NULL) +
    ggplot2::scale_size_area(max_size = 4, limits = c(0, cap), name = "|residual|")
  if(means && all(c("pobs", "pexp", "n") %in% names(d0))){
    # mean length at age of the aged fish, observed and expected under the same sampling by length bin
    mn <- d0[, .(obs = sum(L*n*pobs)/sum(n*pobs), exp = sum(L*n*pexp)/sum(n*pexp), nf = sum(n*pobs)),
             by = .(year, unit, age)][is.finite(obs) & is.finite(exp) & nf >= 1]
    mn[, xa := age + off(unit)]
    p <- p + ggplot2::geom_segment(data = mn, ggplot2::aes(x = xa - 0.18, xend = xa + 0.18, y = exp, yend = exp),
                                   colour = "black", linewidth = 0.6, inherit.aes = FALSE) +
      ggplot2::geom_point(data = mn, ggplot2::aes(xa, obs), shape = 23, fill = "gold", colour = "black", size = 1.6,
                          stroke = 0.4, inherit.aes = FALSE)
  }
  cap_txt <- paste(c(if(!is.null(lhpar) && band) "band: +/- 2 SD of length at age",
                     if(means) "diamonds: observed, bars: expected mean length at age of the aged fish"), collapse = "; ")
  p + ggplot2::facet_wrap(~ year, ncol = ncol) +
    ggplot2::labs(x = "Age", y = "Length", title = fl[1], caption = if(nzchar(cap_txt)) cap_txt,
                  subtitle = paste(attr(x, "type") %||% "osa", "residuals, conditional age-at-length")) +
    theme_ss3() + ggplot2::theme(legend.title = ggplot2::element_text(size = 8))
  })(), text)
}
