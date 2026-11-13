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
#' @param quant \code{"len"} (default) or \code{"age"}
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
ss3_resid <- function(ss3rep, quant = c("len", "age"), type = c("osa", "pearson"), fleets = NULL, seed = 1,
                      long = FALSE){
  quant <- match.arg(quant); type <- match.arg(type)
  db <- as.data.table(if(quant == "len") ss3rep$lendbase else ss3rep$agedbase)
  if(!nrow(db)) stop("no ", quant, " compositions in this SS3 report", call. = FALSE)
  fleets <- .pick_fleets(ss3rep, fleets, .fleetnames(ss3rep)[sort(unique(db$Fleet))])
  db <- db[.fleetnames(ss3rep)[Fleet] %in% fleets]
  ns <- ss3rep$nseasons %||% 1; na <- ss3rep$nareas %||% 1; nsx <- ss3rep$nsexes %||% 1
  dm <- .col(db, c("Nsamp_DM", "DM_effN"))
  d <- data.table(fleet = db$Fleet, name = .fleetnames(ss3rep)[db$Fleet], q = db$Bin, year = db$Yr,
                  unit = .unit_lab(.col(db, "Sex", 1), .col(db, "Sexes", 0), nsx), sex = .col(db, "Sex", 1),
                  season = .season_lab(db$Seas, ns), area = .area_lab(.col(db, "Area", 1), na), iter = "1",
                  part = .col(db, "Part", 0), month = .col(db, "Month", 0), sexes = .col(db, "Sexes", 0),
                  ageerr = .col(db, "Ageerr", 0), repl = .col(db, "Repl.", 1),
                  pobs = db$Obs, pexp = db$Exp, pearson = .col(db, "Pearson", NA_real_), N = db$Nsamp_adj,
                  neff = dm %||% NA_real_)
  d[, partname := name]
  d[, np := uniqueN(part), by = name][np > 1, partname := paste0(name, ".", .part_lab[as.character(part)])][, np := NULL]
  # one composition = fleet x year x season x month x partition x sex option x ageing error x replicate
  d[, comp := .GRP, by = .(fleet, year, season, area, month, part, sexes, ageerr, repl)]
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
  setnames(d, "q", quant)
  # several compositions in one FLQuant cell (months within a season, replicates): Stouffer combination
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
