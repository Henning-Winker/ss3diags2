# resid.R - Pearson and one-step-ahead (OSA) residuals of the SS3 compositions

# beta-binomial cdf at x for N trials, shapes a, b (as FLicc)
.pbetabinom <- function(x, N, a, b){
  if(N <= 0) return(rep(1, length(x)))
  k <- 0:N
  lp <- lchoose(N, k) + lbeta(k + a, N - k + b) - lbeta(a, b)
  cdf <- cumsum(exp(lp - max(lp))); cdf <- cdf/cdf[length(cdf)]
  ifelse(x < 0, 0, ifelse(x >= N, 1, cdf[pmin(floor(x), N) + 1]))
}

# OSA residuals of one composition: counts o, expected proportions p, DM precision a0 (Inf = multinomial)
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
#' @param seed random seed for the OSA randomisation
#' @param long return the long data.table instead of FLQuants
#' @return FLQuants of residuals (one per fleet and partition, as \code{\link{ss3_lencomp}()}: len/age x year x unit x
#'   season x area), or a data.table with the FLQuant dimension columns, \code{data} (residual), \code{obs}, \code{exp}
#'   (counts), \code{n}, \code{a0}, \code{fleet}, \code{name}, \code{partname}; attribute \code{type}
#' @examples
#' data(ss3simple)
#' r <- ss3_resid(ss3simple)
#' names(r)
#' @export
ss3_resid <- function(ss3rep, quant = c("len", "age", "con", "size"), type = c("osa", "pearson"), fleets = NULL,
                      seed = 1, long = FALSE){
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
                  pobs = db$Obs, pexp = db$Exp, pearson = .col(db, "Pearson", NA_real_), N = db$Nsamp_adj,
                  neff = dm %||% NA_real_)
  d[, partname := name]
  d[, np := uniqueN(part), by = name][np > 1, partname := paste0(name, ".", .part_lab[as.character(part)])][, np := NULL]
  # one composition = fleet x year x season x month x partition x sex option x ageing error x replicate
  # (conditional age-at-length: one composition per length bin; size: per size method)
  d[, comp := .GRP, by = .(fleet, year, season, area, month, part, sexes, ageerr, repl, lbin, method)]
  setorder(d, comp, sex, q)
  d[, a0 := { ne <- neff[1]; n <- N[1]
              if(is.finite(ne) && ne > 1 && ne < n) n*(ne - 1)/(n - ne) else Inf }, by = comp]
  d[, `:=`(obs = round(pobs*N), exp = pexp*N, n = N)]
  if(type == "osa"){
    old <- if(exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
    on.exit(if(!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
    d[, data := .osa_comp(obs, pexp, a0[1]), by = comp]
  } else d[, data := pearson]
  setnames(d, "q", qn)
  if(quant == "con" && long){
    d <- d[, c("fleet", "name", "partname", "lbin", qn, .flqdims, "data", "obs", "exp", "n", "a0"), with = FALSE]
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
  out <- long2flqs(d[is.finite(data)], quant = quant, units = paste(type, "residual"), by = "partname",
                   order = unique(d$partname))
  attr(out, "type") <- type
  out
}

utils::globalVariables(c("comp", "a0", "neff", "N", "pobs", "repl", "np"))

#' Plot residuals of conditional age-at-length data
#'
#' OSA (or Pearson) residuals of conditional age-at-length compositions by length bin (x) and age (y), one panel per
#' year, for one fleet: bubbles sized by the absolute residual, blue for observed > expected, red for observed <
#' expected, as the r4ss conditional age-at-length plots. Growth curves (\code{lhpar}) show the expected age at length.
#'
#' @param x long data.table from \code{ss3_resid(ss3rep, quant = "con", long = TRUE)}
#' @param fleet fleet name (default: the first fleet with data)
#' @param years years to plot (default: all with data)
#' @param lhpar optional FLPar from \code{\link{ss3_lhpar}()} for the growth curves by unit
#' @param cap residuals capped at this absolute value for the bubble size
#' @param ncol number of panel columns
#' @return ggplot
#' @export
ss3plot_caal <- function(x, fleet = NULL, years = NULL, lhpar = NULL, cap = 4, ncol = NULL){
  if(!is.data.frame(x) || !"lbin" %in% names(x))
    stop("x must be the long table of ss3_resid(ss3rep, quant = \"con\", long = TRUE)", call. = FALSE)
  d <- as.data.table(x)[is.finite(data)]
  fl <- fleet %||% d$partname[1]
  d <- d[partname %in% fl | name %in% fl]
  if(!is.null(years)) d <- d[year %in% years]
  if(!nrow(d)) stop("no residuals for this fleet / years", call. = FALSE)
  d[, `:=`(age = as.numeric(as.character(age)), r = pmin(abs(data), cap), sign = ifelse(data > 0, "obs > exp", "obs < exp"))]
  p <- ggplot2::ggplot(d, ggplot2::aes(lbin, age)) +
    ggplot2::geom_point(ggplot2::aes(size = r, fill = sign), shape = 21, colour = "grey30", stroke = 0.2, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = c(`obs > exp` = "#3C78B4", `obs < exp` = "#D2604B"), name = NULL) +
    ggplot2::scale_size_area(max_size = 4, limits = c(0, cap), name = "|residual|")
  if(!is.null(lhpar)){
    gc <- .growth_curves(lhpar, maxage = max(d$age))
    gc <- gc[L >= min(d$lbin) & L <= max(d$lbin)*1.05]
    p <- p + ggplot2::geom_line(data = gc, ggplot2::aes(L, age, linetype = unit), colour = "grey40", inherit.aes = FALSE)
  }
  p + ggplot2::facet_wrap(~ year, ncol = ncol) +
    ggplot2::labs(x = "Length bin", y = "Age", title = fl[1],
                  subtitle = paste(attr(x, "type") %||% "osa", "residuals, conditional age-at-length")) +
    theme_ss3() + ggplot2::theme(legend.title = ggplot2::element_text(size = 8))
}
