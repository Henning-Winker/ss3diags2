# plots_growth.R - growth curves by sex and ridge plots of compositions with cohort growth tracks (as FLicc)

#' Mean length and SD at age by unit
#'
#' Mean length and SD at age by unit: from SS3 endgrowth (exact) or from an FLPar (VB with CV linear in length)
#'
#' @param x SS3 report (endgrowth, exact) or FLPar (Schnute or von Bertalanffy, CV linear in length)
#' @param maxage maximum age
#' @param gp growth pattern of an SS3 report (default: the one with the largest share of the recruitment)
#' @return data.table unit, age, L, sd
#' @noRd
.growth_curves <- function(x, maxage = NULL, gp = NULL){
  if(.is_ss3rep(x)){
    # all seasons: mean length at the middle of each season (one point per year with one season), for the morphs
    # of the growth pattern (and its settlement), as ss3_lhpar()
    eg <- as.data.table(x$endgrowth)
    mo <- .gp_morph(x, gp)
    if(!is.null(mo$morph) && "Morph" %in% names(eg) && any(eg$Morph %in% mo$morph)) eg <- eg[Morph %in% mo$morph]
    eg <- eg[!duplicated(eg[, .(Sex, Age_Mid)])][order(Sex, Age_Mid)]
    eg <- eg[Age_Mid >= mo$settle]                                         # not yet settled: no length
    nsx <- x$nsexes %||% 1
    return(eg[, .(unit = if(nsx == 1) "unique" else c("F", "M")[Sex], age = Age_Mid, L = Len_Mid, sd = SD_Mid)])
  }
  units <- if("unit" %in% names(dimnames(x))) dimnames(x)$unit else "unique"
  rbindlist(lapply(units, function(u){
    vb <- .vb_unit(x, u)
    g <- function(p) if(p %in% dimnames(x)$params) { a <- x@.Data; if(length(dim(a)) == 3) a[p, if(u %in% units) u else 1, 1] else a[p, 1] } else NA
    am <- maxage %||% max(vb[["t0"]] + log(20)/vb[["k"]], g("amax") + 1, na.rm = TRUE)
    age <- seq(0, am, length.out = 200)
    cv1 <- g("CVmin"); cv2 <- g("CVmax"); L1 <- g("L1"); L2 <- g("L2")
    A1 <- g("A1"); A2 <- g("A2"); K <- g("K"); bb <- g("richards")
    st <- g("settle"); if(!is.finite(st)) st <- 0
    L <- if(all(is.finite(c(A1, A2, L1, L2, K, bb))) && L1 > 0){
      # SS3 growth, Schnute (1981): Richards (bb != 1) or von Bertalanffy (bb = 1) between L1 at A1 and L2 at A2,
      # linear from 0 below A1; A1 and A2 are ages since settlement, the curve is drawn against calendar age
      ga <- age - st
      fr <- (1 - exp(-K*(ga - A1)))/(1 - exp(-K*(A2 - A1)))
      s <- suppressWarnings(if(abs(bb) > 1e-8) (L1^bb + (L2^bb - L1^bb)*fr)^(1/bb) else L1*exp(log(L2/L1)*fr))
      ifelse(ga < A1, pmax(L1*ga/A1, 0), s)
    } else pmax(vb[["linf"]]*(1 - exp(-vb[["k"]]*(age - vb[["t0"]]))), 0)
    # CV (or SD) of length at age, linear between A1 and A2 in length (f(LAA)) or age (f(A)), as SS3; A2 = 999
    # stands for the accumulator age
    ct <- g("cvtype"); if(!is.finite(ct)) ct <- 0
    if(ct %in% c(1, 3) && all(is.finite(c(A1, A2)))){
      A2c <- if(A2 >= 999) g("amax") else A2
      fr <- pmin(pmax((age - st - A1)/(A2c - A1), 0), 1)
    } else fr <- if(all(is.finite(c(L1, L2)))) pmin(pmax((L - L1)/(L2 - L1), 0), 1) else NA_real_
    v <- if(all(is.finite(c(cv1, cv2)))) cv1 + (cv2 - cv1)*fr else NA_real_
    sa <- g("sdadd"); if(!is.finite(sa)) sa <- 0
    data.table(unit = u, age = age, L = L, sd = (if(ct >= 2) v else v*L) + sa)
  }))
}

#' Plot growth curves by sex
#'
#' Mean length at age with \eqn{\pm} 2 SD bands by sex (unit; r4ss colours, F red, M blue, combined green; a
#' two-sex model with the same growth for both sexes is drawn as one grey curve, "F = M"), with
#' Linf (dashed) and length at 50\% maturity (dotted), and optionally the largest observed length by fleet, as
#' \code{FLicc::plot_growth()}. From an SS3 report the SS3 mean length and SD at age (\code{endgrowth}) are used; from
#' an \code{FLPar} (\code{\link{ss3_lhpar}()}, or FLicc \code{lhpar}) the von Bertalanffy curve with the CV of length
#' at age interpolated between \code{CVmin} (at L1) and \code{CVmax} (at L2).
#'
#' @param x SS3 report or FLPar (\code{linf}, \code{k}, \code{t0}; optional \code{CVmin}, \code{CVmax}, \code{L1},
#'   \code{L2}, \code{L50}; with or without a unit dimension)
#' @param lhpar optional FLPar for Linf and L50 lines when \code{x} is an SS3 report (default \code{ss3_lhpar(x)})
#' @param lfd optional FLQuants of length compositions: largest observed length by fleet as reference lines
#' @param maxage maximum age shown
#' @param gp growth pattern (models with several growth patterns; default: the one with the largest share of the
#'   recruitment)
#' @param text show the title, subtitle and caption (\code{FALSE}: none; labels inside the plot and legends stay)
#' @return ggplot
#' @examples
#' data(ss3simple)
#' ss3plot_growth(ss3simple, lfd = ss3_lencomp(ss3simple))
#' ss3plot_growth(ss3_lhpar(ss3simple))
#' hke <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18.rds", package = "ss3diags2"))
#' ss3plot_growth(hke)                               # two sexes, Richards (Schnute) growth
#' @export
ss3plot_growth <- function(x, lhpar = NULL, lfd = NULL, maxage = NULL, gp = NULL, text = TRUE){
  .drop_text((function(){
  d <- .growth_curves(x, maxage, gp)
  lh <- lhpar %||% (if(.is_ss3rep(x)) suppressMessages(ss3_lhpar(x, gp = gp)) else x)
  if(!is.null(maxage)) d <- d[age <= maxage]
  # two-sex model with the same growth for both sexes (e.g. male offsets of zero): one curve, labelled "F = M"
  same <- FALSE
  if(all(c("F", "M") %in% d$unit)){
    a <- d[unit == "F"]; b <- d[unit == "M"]
    same <- nrow(a) == nrow(b) && isTRUE(all.equal(a$L, b$L, tolerance = 1e-6)) &&
      isTRUE(all.equal(a$sd, b$sd, tolerance = 1e-6))
    if(same) d <- d[unit == "F"]
  }
  units <- unique(d$unit)
  ref <- rbindlist(lapply(units, function(u){
    vb <- .vb_unit(lh, u)
    l50 <- tryCatch({ a <- lh@.Data; if(length(dim(a)) == 3) a["L50", 1, 1] else a["L50", 1] }, error = function(e) NA)
    data.table(unit = u, what = c("Linf", "L50"), L = c(vb[["linf"]], if(u %in% c("F", "unique")) l50 else NA))
  }))[is.finite(L)]
  cols <- ss3unitcols[units]
  if(same){
    d[, unit := "F = M"]; ref[, unit := "F = M"]
    cols <- c(`F = M` = "grey30")
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(age, L, colour = unit, fill = unit))
  if(any(is.finite(d$sd))) p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(L - 2*sd, 0), ymax = L + 2*sd),
                                                        colour = NA, alpha = 0.15)
  p <- p + ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(data = ref, ggplot2::aes(yintercept = L, colour = unit, linetype = what), linewidth = 0.5) +
    ggplot2::scale_linetype_manual(values = c(Linf = 2, L50 = 3, `largest fish` = 1), name = NULL)
  if(!is.null(lfd)){
    ml <- .flq_long(lfd, 1, "lfd")[is.finite(data) & data > 0, .(L = max(as.numeric(as.character(quant)))), by = name]
    ml <- ml[, .(name = paste(unique(name), collapse = ", ")), by = L]        # one label for fleets with the same max
    p <- p + ggplot2::geom_hline(data = ml, ggplot2::aes(yintercept = L), colour = "grey50", linewidth = 0.3) +
      ggplot2::annotate("text", x = 0, y = ml$L, label = paste("max", ml$name), hjust = 0, vjust = -0.3, size = 2.5,
                        colour = "grey40")
  }
  p + ggplot2::scale_colour_manual(values = cols, name = NULL) + ggplot2::scale_fill_manual(values = cols, name = NULL) +
    ggplot2::expand_limits(y = 0) + ggplot2::labs(x = "Age (years)", y = "Length") + theme_ss3()
  })(), text)
}

#' Ridge plots of length or age compositions
#'
#' Compositions as ridges along the year axis (one ridge per year, proportions scaled within each panel), one panel per
#' fleet and unit (sex), season or area with data, as \code{FLicc::plot_growth_lfd()}: observed (filled), expected
#' (filled) or both (observed filled, expected as lines), with cohort tracks of von Bertalanffy growth by sex from
#' \code{lhpar} (length; diagonals for ages) to follow year classes through the compositions.
#'
#' @param comps FLQuants of observed numbers (or proportions) at length or age (e.g. \code{\link{ss3_lencomp}()})
#' @param fit FLQuants of expected compositions on the same scale (e.g. \code{ss3_lencomp(ss3rep, what = "exp")})
#' @param type \code{"obs"} (default), \code{"exp"} or \code{"both"} (expected as lines over the observed)
#' @param by dimensions kept as separate panels, any of \code{"unit"} (default), \code{"season"}, \code{"area"}
#' @param lhpar FLPar for the cohort tracks (length compositions; \code{NULL} for none); ages get diagonals
#' @param cohort_by years between cohort tracks
#' @param years optional years
#' @param scale ridge height relative to the year spacing
#' @param ncol facet columns
#' @param fleets fleets to plot: names (a fleet name selects all its partitions, e.g. \code{"DK3a"} gives
#'   \code{DK3a.discard} and \code{DK3a.retained}) or positions in the FLQuants (default: all)
#' @param units units (sexes) to plot, e.g. \code{"F"}, \code{"M"} or \code{"unique"} (combined-sex data); fleets
#'   without these units are dropped (default: all)
#' @param text show the title, subtitle and caption (\code{FALSE}: none; labels inside the plot and legends stay)
#' @return ggplot
#' @examples
#' data(ss3simple)
#' lfd <- ss3_lencomp(ss3simple)
#' ss3plot_ridges(lfd["SURVEY1"], fit = ss3_lencomp(ss3simple, what = "exp")["SURVEY1"], type = "both",
#'                lhpar = ss3_lhpar(ss3simple))
#' @export
ss3plot_ridges <- function(comps, fit = NULL, type = c("obs", "exp", "both"), by = "unit", lhpar = NULL,
                           cohort_by = 2, years = NULL, scale = 0.9, ncol = NULL,
                           fleets = NULL, units = NULL, text = TRUE){
  .drop_text((function(){
  comps <- .subset_flqs(comps, fleets, units)
  if(!is.null(fit)) fit <- .subset_flqs(fit, fleets, units)
  type <- match.arg(type, c("obs", "exp", "both")); dims <- .bydims(by)
  if(type != "obs" && is.null(fit)) stop("type = \"", type, "\" needs 'fit'", call. = FALSE)
  qn <- FLCore::quant(if(is(comps, "FLQuant")) comps else comps[[1]])
  prep <- function(x){
    d <- .collapse(.flq_long(x, 1, "comps"), dims)[is.finite(data)]
    d[, q := as.numeric(as.character(quant))]
    if(!is.null(years)) d <- d[year %in% years]
    d[, tot := sum(data), by = .(name, year, unit, season, area)][tot > 0]
    d[, p := data/tot]
    d <- .add_group(d, dims); .panel(d, dims)
  }
  o <- prep(comps); e <- if(!is.null(fit)) prep(fit) else NULL
  if(!is.null(e)){ e <- e[paste(panel, year) %in% paste(o$panel, o$year)]; e[, panel := factor(as.character(panel), levels(o$panel))] }
  bw <- .bw(o$q)/0.9
  # each panel on the observed range (first to last bin with observations)
  rg <- o[p > 0, .(qlo = min(q), qhi = max(q)), by = panel]
  trim <- function(d) { d <- merge(d, rg, by = "panel"); d[q >= qlo & q <= qhi][, c("qlo", "qhi") := NULL] }
  o <- trim(o); if(!is.null(e)) e <- trim(e)
  mx <- rbind(o, if(!is.null(e)) e)[, .(mx = max(p)), by = panel]
  ridge <- function(d){
    d <- merge(d, mx, by = "panel")[order(panel, year, q)]
    d[, `:=`(x = year + scale*p/mx, y = q + if(qn == "len") bw/2 else 0)]
    d
  }
  poly <- function(d) rbindlist(lapply(split(d, by = c("panel", "year"), drop = TRUE), function(z)
    rbind(z[, .(panel, year, x, y)], z[.N:1, .(panel, year, x = year, y)])))
  ob <- ridge(o); ex <- if(!is.null(e)) ridge(e) else NULL
  fillc <- "#6BAED6"
  p <- ggplot2::ggplot()
  if(type %in% c("obs", "both")) p <- p + ggplot2::geom_polygon(data = poly(ob), ggplot2::aes(x, y, group = interaction(panel, year)),
                                                                fill = fillc, colour = NA, alpha = 0.8)
  if(type == "exp") p <- p + ggplot2::geom_polygon(data = poly(ex), ggplot2::aes(x, y, group = interaction(panel, year)),
                                                   fill = "grey60", colour = NA, alpha = 0.8)
  if(type == "both") p <- p + ggplot2::geom_path(data = ex, ggplot2::aes(x, y, group = interaction(panel, year)),
                                                 colour = "black", linewidth = 0.35)
  if(qn == "age" || !is.null(lhpar)){
    tr <- .cohort_lines(o[, .(panel, unit, year, q)], qn, lhpar, cohort_by, binw = if(qn == "len") bw else 0)
    if(nrow(tr)) p <- p + ggplot2::geom_line(data = tr, ggplot2::aes(year, q, group = interaction(birth, lunit), linetype = lunit),
                                             colour = "#B2182B", linewidth = 0.3) +
      ggplot2::scale_linetype_manual(values = c(F = 2, M = 3, unique = 2),
                                     guide = if(length(unique(tr$lunit)) > 1) ggplot2::guide_legend(title = "growth") else "none")
  }
  p + ggplot2::facet_wrap(~ panel, scales = "free", ncol = ncol) +
    ggplot2::labs(x = "Year", y = switch(qn, len = "Length", age = "Age", qn),
                  subtitle = switch(type, obs = "observed", exp = "expected", both = "observed (filled), expected (lines)")) +
    scale_x_year() + theme_ss3()
  })(), text)
}

utils::globalVariables(c("qlo", "qhi", "tot", "p", "mx", "y", "what", "Len_Mid", "SD_Mid", "Age_Mid", "L", "sd"))

#' Plot sex ratios of compositions
#'
#' Observed (points) and expected (line) proportion of females by length (or age) bin, pooled over years
#' (\code{type = "len"}), or by year, pooled over bins (\code{type = "year"}), for the fleets with sex-specific
#' compositions (units "F" and "M"). The proportions use the numbers of the compositions (\code{scale = "n"}), so
#' years and bins are weighted by their sample sizes.
#'
#' @param comps FLQuants of compositions with units F and M (\code{ss3_lencomp()}, \code{ss3_agecomp()})
#' @param fit optional FLQuants of expected compositions (\code{what = "exp"})
#' @param type \code{"len"} (by bin, default) or \code{"year"}
#' @param fleets,units fleets to plot (as \code{\link{ss3plot_comps}()}; units are F and M)
#' @param ncol number of panel columns
#' @param ci level of the binomial confidence interval of the observed proportions (Wilson score interval with the
#'   numbers in the bin, i.e. the input or adjusted sample size of the compositions); \code{NULL} or \code{FALSE}
#'   for none. The interval reflects sampling error of the observations; the SS3 hessian gives no standard errors for
#'   the expected compositions
#' @param min.n minimum (effective) number of fish in a bin for an interval; bins with fewer, mostly the largest
#'   lengths, are nearly uninformative and their intervals span most of 0-1
#' @param text show the title, subtitle and caption (\code{FALSE}: none; labels inside the plot and legends stay)
#' @return ggplot
#' @examples
#' hke <- readRDS(system.file("extdata", "hke.37.17-18", "hke.37.17-18.rds", package = "ss3diags2"))
#' ss3plot_sexratio(ss3_lencomp(hke), fit = ss3_lencomp(hke, what = "exp"))
#' ss3plot_sexratio(ss3_lencomp(hke), fit = ss3_lencomp(hke, what = "exp"), type = "year", ci = FALSE)
#' @export
ss3plot_sexratio <- function(comps, fit = NULL, type = c("len", "year"), fleets = NULL, ncol = NULL, ci = 0.95, min.n = 5, text = TRUE){
  .drop_text((function(){
  type <- match.arg(type, c("len", "year"))
  comps <- .subset_flqs(comps, fleets, NULL)
  keep <- names(comps)[vapply(comps, function(q) all(c("F", "M") %in% dimnames(q)$unit), TRUE)]
  if(!length(keep)) stop("no compositions with units F and M", call. = FALSE)
  sr <- function(x, what){
    d <- .flq_long(x[intersect(keep, names(x))], 1, what)[unit %in% c("F", "M") & is.finite(data)]
    d[, q := as.numeric(as.character(quant))]
    by <- c("name", if(type == "len") "q" else "year")
    d <- dcast(d[, .(v = sum(data)), by = c(by, "unit")], stats::as.formula(paste(paste(by, collapse = "+"), "~ unit")),
               value.var = "v")
    d[, `:=`(pF = F/(F + M), n = F + M, what = what)][is.finite(pF)]
  }
  o <- sr(comps, "observed")
  p <- ggplot2::ggplot(o, ggplot2::aes(if(type == "len") q else year, pF)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = 3, colour = "grey50")
  if(!is.null(ci) && !isFALSE(ci) && ci > 0){
    # Wilson score interval for a binomial proportion with the (effective) sample size n in the bin / year
    z <- stats::qnorm(1 - (1 - ci)/2)
    o[, `:=`(lo = (pF + z^2/(2*n) - z*sqrt(pF*(1 - pF)/n + z^2/(4*n^2)))/(1 + z^2/n),
             hi = (pF + z^2/(2*n) + z*sqrt(pF*(1 - pF)/n + z^2/(4*n^2)))/(1 + z^2/n))]
    p <- p + ggplot2::geom_linerange(data = o[n >= min.n], ggplot2::aes(ymin = lo, ymax = hi), colour = "grey55",
                                     linewidth = 0.35)
  }
  p <- p +
    ggplot2::geom_point(ggplot2::aes(size = n), shape = 21, fill = "grey60", colour = "grey20", alpha = 0.8) +
    ggplot2::scale_size_area(max_size = 3, guide = "none")
  if(!is.null(fit)){
    e <- sr(fit, "expected")
    # expected composition scaled like the observed: use the observed totals by bin / year as weights
    p <- p + ggplot2::geom_line(data = e, ggplot2::aes(if(type == "len") q else year, pF), colour = ss3unitcols[["F"]],
                                linewidth = 0.8, inherit.aes = FALSE)
  }
  p + ggplot2::facet_wrap(~ name, ncol = ncol) + ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::labs(x = if(type == "len") "Length" else "Year", y = "Proportion female",
                  caption = paste0("points: observed (size: numbers)", if(!is.null(ci) && !isFALSE(ci) && ci > 0)
                    paste0(", bars: ", round(100*ci), "% binomial CI"), if(!is.null(fit)) "; line: expected")) +
    theme_ss3()
  })(), text)
}


#' Growth pattern, settlement and morphs of an SS3 model
#'
#' @param ss3rep SS3 report (\code{r4ss::SS_output()})
#' @param gp growth pattern (default: the one with the largest share of the recruitment)
#' @return list: gp, settle (time of settlement as a fraction of the year), sseas (settlement season) and morph
#'   (morph index of females and males, or \code{NULL} without morph indexing)
#' @noRd
.gp_morph <- function(ss3rep, gp = NULL){
  rd <- ss3rep$recruitment_dist$recruit_dist
  settle <- 0; sseas <- 1L
  if(is.data.frame(rd) && "Settle_Month" %in% names(rd) && nrow(rd)){
    w <- rowSums(rd[, grep("^recr_dist", names(rd)), drop = FALSE], na.rm = TRUE)
    rgp <- if("G_pattern" %in% names(rd)) as.integer(as.character(rd$G_pattern)) else rep(1L, nrow(rd))
    if(is.null(gp)) gp <- rgp[if(any(w > 0)) which.max(w) else 1]
    i <- which(rgp == gp); if(!length(i)) i <- seq_along(rgp)
    i <- i[which.max(w[i])]
    settle <- (suppressWarnings(as.numeric(as.character(rd$Settle_Month[i]))) - 1)/12
    sseas <- suppressWarnings(as.integer(as.character(rd$Seas[i])))
    if(length(settle) != 1 || !is.finite(settle) || settle < 0) settle <- 0
  }
  gp <- gp %||% 1L
  mi <- ss3rep$morph_indexing; morph <- NULL
  if(is.data.frame(mi) && all(c("Index", "GP", "Sex", "BirthSeas") %in% names(mi)))
    morph <- mi$Index[mi$GP == gp & mi$BirthSeas == sseas]
  list(gp = gp, settle = settle, sseas = sseas, morph = morph)
}
