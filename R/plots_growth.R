# plots_growth.R - growth curves by sex and ridge plots of compositions with cohort growth tracks (as FLicc)

# mean length and SD at age by unit: from SS3 endgrowth (exact) or from an FLPar (VB with CV linear in length)
.growth_curves <- function(x, maxage = NULL){
  if(.is_ss3rep(x)){
    eg <- as.data.table(x$endgrowth)[Seas == 1]
    eg <- eg[!duplicated(eg[, .(Sex, Age_Beg)])]
    nsx <- x$nsexes %||% 1
    return(eg[, .(unit = if(nsx == 1) "unique" else c("F", "M")[Sex], age = Age_Mid, L = Len_Mid, sd = SD_Mid)])
  }
  units <- if("unit" %in% names(dimnames(x))) dimnames(x)$unit else "unique"
  rbindlist(lapply(units, function(u){
    vb <- .vb_unit(x, u)
    g <- function(p) if(p %in% dimnames(x)$params) { a <- x@.Data; if(length(dim(a)) == 3) a[p, if(u %in% units) u else 1, 1] else a[p, 1] } else NA
    am <- maxage %||% (vb[["t0"]] + log(20)/vb[["k"]])
    age <- seq(0, am, length.out = 200)
    L <- pmax(vb[["linf"]]*(1 - exp(-vb[["k"]]*(age - vb[["t0"]]))), 0)
    cv1 <- g("CVmin"); cv2 <- g("CVmax"); L1 <- g("L1"); L2 <- g("L2")
    cv <- if(all(is.finite(c(cv1, cv2, L1, L2)))) cv1 + (cv2 - cv1)*pmin(pmax((L - L1)/(L2 - L1), 0), 1) else NA_real_
    data.table(unit = u, age = age, L = L, sd = cv*L)
  }))
}

#' Plot growth curves by sex
#'
#' Mean length at age with \eqn{\pm} 2 SD bands by sex (unit; r4ss colours, F red, M blue, combined green), with
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
#' @return ggplot
#' @examples
#' data(ss3ank)
#' ss3plot_growth(ss3ank, lfd = ss3_lencomp(ss3ank))
#' ss3plot_growth(ss3_lhpar(ss3ank))
#' @export
ss3plot_growth <- function(x, lhpar = NULL, lfd = NULL, maxage = NULL){
  d <- .growth_curves(x, maxage)
  lh <- lhpar %||% (if(.is_ss3rep(x)) ss3_lhpar(x) else x)
  if(!is.null(maxage)) d <- d[age <= maxage]
  units <- unique(d$unit)
  ref <- rbindlist(lapply(units, function(u){
    vb <- .vb_unit(lh, u)
    l50 <- tryCatch({ a <- lh@.Data; if(length(dim(a)) == 3) a["L50", 1, 1] else a["L50", 1] }, error = function(e) NA)
    data.table(unit = u, what = c("Linf", "L50"), L = c(vb[["linf"]], if(u %in% c("F", "unique")) l50 else NA))
  }))[is.finite(L)]
  cols <- ss3unitcols[units]
  p <- ggplot2::ggplot(d, ggplot2::aes(age, L, colour = unit, fill = unit))
  if(any(is.finite(d$sd))) p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = pmax(L - 2*sd, 0), ymax = L + 2*sd),
                                                        colour = NA, alpha = 0.15)
  p <- p + ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(data = ref, ggplot2::aes(yintercept = L, colour = unit, linetype = what), linewidth = 0.5) +
    ggplot2::scale_linetype_manual(values = c(Linf = 2, L50 = 3, `largest fish` = 1), name = NULL)
  if(!is.null(lfd)){
    ml <- .flq_long(lfd, 1, "lfd")[is.finite(data) & data > 0, .(L = max(as.numeric(as.character(quant)))), by = name]
    p <- p + ggplot2::geom_hline(data = ml, ggplot2::aes(yintercept = L), colour = "grey50", linewidth = 0.3) +
      ggplot2::annotate("text", x = 0, y = ml$L, label = paste("max", ml$name), hjust = 0, vjust = -0.3, size = 2.5,
                        colour = "grey40")
  }
  p + ggplot2::scale_colour_manual(values = cols, name = NULL) + ggplot2::scale_fill_manual(values = cols, name = NULL) +
    ggplot2::expand_limits(y = 0) + ggplot2::labs(x = "Age (years)", y = "Length") + theme_ss3()
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
#' @return ggplot
#' @examples
#' data(ss3ank)
#' lfd <- ss3_lencomp(ss3ank)["FR_IE_IBTS"]
#' ss3plot_ridges(lfd, fit = ss3_lencomp(ss3ank, what = "exp")["FR_IE_IBTS"], type = "both",
#'                lhpar = ss3_lhpar(ss3ank))
#' @export
ss3plot_ridges <- function(comps, fit = NULL, type = c("obs", "exp", "both"), by = "unit", lhpar = NULL,
                           cohort_by = 2, years = NULL, scale = 0.9, ncol = NULL,
                           fleets = NULL, units = NULL){
  comps <- .subset_flqs(comps, fleets, units)
  if(!is.null(fit)) fit <- .subset_flqs(fit, fleets, units)
  type <- match.arg(type); dims <- .bydims(by)
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
}

utils::globalVariables(c("qlo", "qhi", "tot", "p", "mx", "y", "what", "Len_Mid", "SD_Mid", "Age_Mid", "L", "sd"))
