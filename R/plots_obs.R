# plots_obs.R - plots of observations (FLR inputs)

#' Plot catch by fleet
#'
#' Catch by fleet from an \code{FLQuants} (one \code{FLQuant} per fleet, e.g. \code{\link{ss3_catch}()} or the
#' JABBA3 / FLicc catch input), as stacked bars, lines or shares, optionally split by unit, season or area.
#'
#' @param catch FLQuants (or FLQuant) of catch by fleet
#' @param by dimensions to keep separate, any of \code{"unit"}, \code{"season"}, \code{"area"} (panels by row);
#'   other dimensions are summed
#' @param type \code{"bar"} (stacked, default), \code{"line"} or \code{"share"} (proportion by fleet)
#' @param iter iteration to plot
#' @param ylab y-axis label (default from the units)
#' @return ggplot
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' ss3plot_catch(ss3_catch(ss3sma))
#' ss3plot_catch(ss3_catch(ss3sma), type = "share")
#' }
#' @export
ss3plot_catch <- function(catch, by = NULL, type = c("bar", "line", "share"), iter = 1, ylab = NULL){
  type <- match.arg(type); dims <- .bydims(by)
  d <- .flq_long(catch, iter, "catch")
  un <- unique(vapply(if(is(catch, "FLQuant")) list(catch) else catch, FLCore::units, ""))
  d <- .collapse(d, dims)[!is.na(data)]
  d <- .add_group(d, dims)
  cols <- stats::setNames(ss3cols(nlevels(d$name)), levels(d$name))
  ylab <- ylab %||% if(type == "share") "Share of catch" else
    paste0("Catch (", paste(ifelse(un == "1000", "1000s", un), collapse = " | "), ")")
  p <- ggplot2::ggplot(d, ggplot2::aes(year, data))
  p <- if(type == "line")
    p + ggplot2::geom_line(ggplot2::aes(colour = name), linewidth = 0.6) + ggplot2::scale_colour_manual(values = cols)
  else
    p + ggplot2::geom_col(ggplot2::aes(fill = name), position = if(type == "share") "fill" else "stack",
                          width = 0.9, colour = NA) + ggplot2::scale_fill_manual(values = cols)
  if(length(dims)) p <- p + ggplot2::facet_grid(group ~ ., scales = if(type == "share") "fixed" else "free_y")
  if(length(un) > 1 && type != "share") warning("catch in mixed units: ", paste(un, collapse = ", "), call. = FALSE)
  p + ggplot2::labs(x = "Year", y = ylab) + scale_x_year() + theme_ss3()
}

#' Plot abundance indices
#'
#' Observed indices with lognormal confidence intervals from \code{index.var} (log-scale SE), from \code{FLIndices}
#' (e.g. \code{\link{ss3_indices}()} or the JABBA3 input), with an optional fit (expected index) overlaid, one panel
#' per index, colours by unit, season or area (r4ss season colours).
#'
#' @param indices FLIndices (or FLIndexBiomass)
#' @param fit optional FLIndices (e.g. \code{ss3_indices(ss3rep, what = "exp")}) or FLQuants with the expected
#'   indices, matched by name
#' @param by dimensions shown in colour, any of \code{"unit"}, \code{"season"}, \code{"area"}
#' @param ci confidence level of the intervals (0 for none)
#' @param log plot on log scale
#' @param combine all indices in one panel, each scaled to its mean
#' @param time x axis as year + timing within the year (\code{range} startf, and season for seasonal indices)
#' @param ncol number of panel columns
#' @param iter iteration to plot
#' @return ggplot
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' ss3plot_indices(ss3_indices(ss3sma), fit = ss3_indices(ss3sma, what = "exp"))
#' ss3plot_indices(ss3_indices(ss3sma), combine = TRUE)
#' }
#' @export
ss3plot_indices <- function(indices, fit = NULL, by = NULL, ci = 0.95, log = FALSE, combine = FALSE, time = FALSE,
                            ncol = NULL, iter = 1){
  dims <- .bydims(by)
  d <- .fli_long(indices, iter)
  ns <- max(1, length(setdiff(unique(d$season), "all")))
  d[, xval := if(time) year + ifelse(season == "all", timef, (as.numeric(season) - 1)/ns + timef/ns) else year]
  z <- if(ci > 0) stats::qnorm(1 - (1 - ci)/2) else 0
  d[, `:=`(lo = data*exp(-z*se), hi = data*exp(z*se))]
  if(!is.null(fit)){
    f <- if(is(fit, "FLIndices") || is(fit, "FLIndex") || is(fit, "FLIndexBiomass")) .fli_long(fit, iter) else
      .flq_long(fit, iter, "fit")
    setnames(f, "data", "fit"); lv <- levels(d$name)
    d <- merge(d, f[, .(name = as.character(name), quant, year, unit, season, area, fit)],
               by = c("name", "quant", "year", "unit", "season", "area"), all.x = TRUE, sort = FALSE)
    d[, name := factor(name, levels = lv)]
  }
  d <- .add_group(d, dims)
  if(combine){
    d[, mu := mean(data, na.rm = TRUE), by = name]
    for(v in intersect(c("data", "lo", "hi", "fit"), names(d))) d[, (v) := get(v)/mu]
    cols <- stats::setNames(ss3cols(nlevels(d$name)), levels(d$name))
    p <- ggplot2::ggplot(d, ggplot2::aes(xval, data, colour = name)) +
      ggplot2::geom_line(linewidth = 0.4, alpha = 0.6) + ggplot2::geom_point(size = 1.6) +
      ggplot2::scale_colour_manual(values = cols) + ggplot2::labs(x = "Year", y = "Index / mean")
    if(log) p <- p + ggplot2::scale_y_log10()
    return(p + scale_x_year() + theme_ss3())
  }
  cols <- .group_cols(d, dims)
  p <- ggplot2::ggplot(d, ggplot2::aes(xval, data))
  if(ci > 0) p <- p + ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi, colour = group), linewidth = 0.4,
                                              alpha = 0.7, show.legend = FALSE)
  if(!is.null(fit) && "fit" %in% names(d)){
    df <- d[!is.na(fit)][order(name, group, year)][, seg := paste(group, .gap_run(year)), by = .(name, group)]
    p <- p + ggplot2::geom_line(data = df, ggplot2::aes(y = fit, group = seg), colour = "blue", linewidth = 0.7)
  }
  p <- p + ggplot2::geom_point(ggplot2::aes(fill = group), shape = 21, size = 1.8, colour = "black") +
    ggplot2::scale_colour_manual(values = cols) + ggplot2::scale_fill_manual(values = cols) +
    ggplot2::facet_wrap(~ name, scales = "free_y", ncol = ncol) +
    ggplot2::labs(x = "Year", y = "Index")
  p <- if(log) p + ggplot2::scale_y_log10() else p + ggplot2::expand_limits(y = 0)
  p <- p + scale_x_year() + theme_ss3()
  if(!length(dims)) p <- p + ggplot2::theme(legend.position = "none")
  p
}

#' Plot length or age compositions
#'
#' Compositions by fleet from \code{FLQuants} (one \code{FLQuant} per fleet with length bins or ages as the first
#' dimension, e.g. \code{\link{ss3_lencomp}()}, \code{\link{ss3_agecomp}()} or the FLicc \code{lfd}), optionally split
#' by unit (e.g. sex), season or area; dimensions not in \code{by} are summed (numbers).
#'
#' @param comps FLQuants (or FLQuant) of numbers (or proportions) at length or age
#' @param fit optional FLQuants of expected compositions on the same scale (e.g. \code{ss3_lencomp(ss3rep, what =
#'   "exp")}), shown as lines
#' @param by dimensions kept separate, any of \code{"unit"}, \code{"season"}, \code{"area"} (colours, r4ss
#'   conventions: unique green, F red, M blue; r4ss season colours)
#' @param type \code{"pooled"} (summed over years, default), \code{"annual"} (one panel per year), \code{"bubble"}
#'   (year x bin, size = proportion) or \code{"mean"} (mean length/age by year)
#' @param years optional years to show
#' @param ncol number of panel columns (\code{"annual"} with one fleet)
#' @param iter iteration to plot
#' @return ggplot
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' lfd <- ss3_lencomp(ss3sma)
#' ss3plot_comps(lfd, by = "unit")
#' ss3plot_comps(lfd["Fishery_4"], type = "annual", by = "unit", fit = ss3_lencomp(ss3sma, what = "exp")["Fishery_4"])
#' ss3plot_comps(lfd, type = "bubble")
#' ss3plot_comps(lfd, type = "mean", by = "unit")
#' }
#' @export
ss3plot_comps <- function(comps, fit = NULL, by = NULL, type = c("pooled", "annual", "bubble", "mean"),
                          years = NULL, ncol = NULL, iter = 1){
  type <- match.arg(type); dims <- .bydims(by)
  q0 <- if(is(comps, "FLQuant")) comps else comps[[1]]
  qn <- FLCore::quant(q0)
  xlab <- switch(qn, len = "Length", age = "Age", qn)
  d <- .flq_long(comps, iter, "comps")
  d <- .collapse(d, dims)
  if(!is.null(fit)){
    f <- .collapse(.flq_long(fit, iter, "comps"), dims)
    setnames(f, "data", "fit"); lv <- levels(d$name)
    d <- merge(d, f[, .(name = as.character(name), quant, year, unit, season, area, fit)],
               by = c("name", "quant", "year", "unit", "season", "area"), all.x = TRUE, sort = FALSE)
    d[, name := factor(name, levels = lv)]
  } else d[, fit := NA_real_]
  if(!is.null(years)) d <- d[year %in% years]
  d[, quant := as.numeric(as.character(quant))]
  # drop compositions without data
  d[, ntot := sum(data, na.rm = TRUE), by = .(name, year, unit, season, area)]
  d <- d[ntot > 0 & !is.na(ntot)]
  d <- .add_group(d, dims)
  cols <- .group_cols(d, dims)
  if(!length(dims)) cols[] <- "#008000"
  gk <- c("name", "year", "group")

  if(type == "pooled"){
    a <- d[, .(data = sum(data, na.rm = TRUE), fit = sum(fit, na.rm = TRUE)), by = .(name, group, quant)]
    a[, `:=`(prop = data/sum(data), pexp = fit/sum(fit)), by = .(name, group)]
    p <- ggplot2::ggplot(a, ggplot2::aes(quant, prop, colour = group))
    p <- if(!length(dims)) p + ggplot2::geom_col(fill = "grey75", colour = NA, width = .bw(a$quant)) else
      p + ggplot2::geom_line(linewidth = 0.6)
    if(!is.null(fit)) p <- p + ggplot2::geom_line(ggplot2::aes(y = pexp), linewidth = 0.7,
                                                  linetype = if(length(dims)) 2 else 1)
    p <- p + ggplot2::facet_wrap(~ name, scales = "free_y", ncol = ncol) +
      ggplot2::labs(x = xlab, y = "Proportion (pooled over years)")
  }
  if(type == "annual"){
    d[, `:=`(prop = data/sum(data, na.rm = TRUE), pexp = fit/sum(fit, na.rm = TRUE)), by = gk]
    p <- ggplot2::ggplot(d, ggplot2::aes(quant, prop, colour = group))
    p <- if(!length(dims)) p + ggplot2::geom_col(fill = "grey75", colour = NA, width = .bw(d$quant)) else
      p + ggplot2::geom_line(linewidth = 0.5)
    if(!is.null(fit)) p <- p + ggplot2::geom_line(ggplot2::aes(y = pexp), linewidth = 0.6,
                                                  linetype = if(length(dims)) 2 else 1)
    p <- p + if(nlevels(droplevels(d$name)) == 1) ggplot2::facet_wrap(~ year, ncol = ncol, dir = "v") else
      ggplot2::facet_grid(year ~ name)
    p <- p + ggplot2::labs(x = xlab, y = "Proportion") +
      ggplot2::theme(strip.text = ggplot2::element_text(size = 7), axis.text.y = ggplot2::element_text(size = 6))
  }
  if(type == "bubble"){
    d[, prop := data/sum(data, na.rm = TRUE), by = gk]
    d <- .panel(d[prop > 0], dims)                         # panels with data only, years with data only
    p <- ggplot2::ggplot(d, ggplot2::aes(year, quant, size = prop)) +
      ggplot2::geom_point(shape = 21, fill = ss3bubblecols[["obs"]], colour = "white", stroke = 0.15, alpha = 0.6) +
      ggplot2::scale_size_area(max_size = 3.5, name = "Proportion") +
      ggplot2::facet_wrap(~ panel, scales = "free", ncol = ncol) +
      ggplot2::labs(x = "Year", y = xlab) + scale_x_year(n = 4)
    return(p + theme_ss3() + ggplot2::theme(legend.title = ggplot2::element_text(size = 8)))
  }
  if(type == "mean"){
    a <- d[, .(mu = sum(quant*data, na.rm = TRUE)/sum(data, na.rm = TRUE),
               mexp = if(all(is.na(fit))) NA_real_ else sum(quant*fit, na.rm = TRUE)/sum(fit, na.rm = TRUE)), by = gk]
    p <- ggplot2::ggplot(a, ggplot2::aes(year, mu, colour = group)) +
      ggplot2::geom_point(ggplot2::aes(fill = group), shape = 21, colour = "black", size = 1.8)
    a <- a[order(name, group, year)][, seg := paste(group, .gap_run(year)), by = .(name, group)]
    if(!is.null(fit)) p <- p + ggplot2::geom_line(data = a, ggplot2::aes(y = mexp, group = seg), linewidth = 0.6)
    p <- p + ggplot2::scale_fill_manual(values = cols) +
      ggplot2::facet_wrap(~ name, scales = "free_y", ncol = ncol) +
      ggplot2::labs(x = "Year", y = paste("Mean", tolower(xlab))) + scale_x_year()
  }
  p <- p + ggplot2::scale_colour_manual(values = cols) + theme_ss3()
  if(!length(dims)) p <- p + ggplot2::theme(legend.position = "none")
  p
}

# bin width for bars
.bw <- function(x){ u <- sort(unique(x)); if(length(u) > 1) min(diff(u))*0.9 else 0.9 }

#' Plot composition sample sizes
#'
#' Sample size of each composition (sum over bins of numbers at length or age, i.e. the SS3 sample size when the
#' compositions come from \code{\link{ss3_lencomp}()} with \code{scale = "n"}), by fleet and year, optionally split by
#' unit, season or area.
#' @param comps FLQuants (or FLQuant) of numbers at length or age
#' @param by dimensions kept separate, any of \code{"unit"}, \code{"season"}, \code{"area"}
#' @param iter iteration
#' @return ggplot
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' ss3plot_nsamp(ss3_lencomp(ss3sma), by = "unit")
#' }
#' @export
ss3plot_nsamp <- function(comps, by = NULL, iter = 1){
  dims <- .bydims(by)
  d <- .collapse(.flq_long(comps, iter, "comps"), dims)
  a <- d[, .(data = sum(data, na.rm = TRUE)), by = .(name, year, unit, season, area)][data > 0]
  a <- .add_group(a, dims)
  cols <- .group_cols(a, dims)
  p <- ggplot2::ggplot(a, ggplot2::aes(year, data, fill = group)) +
    ggplot2::geom_segment(ggplot2::aes(xend = year, y = 0, yend = data, colour = group), linewidth = 0.3,
                          show.legend = FALSE,
                          position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::geom_point(shape = 21, size = 1.8, position = ggplot2::position_dodge(width = 0.6)) +
    ggplot2::scale_fill_manual(values = cols) + ggplot2::scale_colour_manual(values = cols) +
    ggplot2::facet_wrap(~ name, scales = "free_y") + ggplot2::expand_limits(y = 0) +
    ggplot2::labs(x = "Year", y = "Sample size") + scale_x_year() + theme_ss3()
  if(!length(dims)) p <- p + ggplot2::theme(legend.position = "none")
  p
}

#' Plot data coverage
#'
#' Which data are available by fleet and year (as \code{r4ss::SSplotData()}), from the FLR objects: catch
#' (\code{FLQuants}), indices (\code{FLIndices}), length and age compositions (\code{FLQuants}). Circle size is
#' relative to the fleet's maximum: catch, index precision (1/SE) and composition sample size.
#'
#' @param catch FLQuants of catch by fleet
#' @param indices FLIndices
#' @param lencomp,agecomp FLQuants of numbers at length / age
#' @param by dimensions kept separate, any of \code{"unit"}, \code{"season"}, \code{"area"} (panels by column)
#' @param iter iteration
#' @return ggplot
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' ss3plot_data(catch = ss3_catch(ss3sma), indices = ss3_indices(ss3sma), lencomp = ss3_lencomp(ss3sma))
#' }
#' @export
ss3plot_data <- function(catch = NULL, indices = NULL, lencomp = NULL, agecomp = NULL, by = NULL, iter = 1){
  dims <- .bydims(by)
  keep <- c("name", "year", "unit", "season", "area", "data", "dtype")
  out <- list()
  if(!is.null(catch)){
    d <- .collapse(.flq_long(catch, iter, "catch"), dims)[!is.na(data) & data > 0]
    out$catch <- d[, dtype := "Catch"][, keep, with = FALSE]
  }
  if(!is.null(indices)){
    d <- .fli_long(indices, iter)
    d[, data := ifelse(is.finite(se) & se > 0, 1/se, 1)]
    d <- .collapse(d[, .(name = as.character(name), year, unit, season, area, data)], dims)
    out$index <- d[, dtype := "Index"][, keep, with = FALSE]
  }
  cmp <- function(x, lab){
    d <- .collapse(.flq_long(x, iter, lab), dims)
    d <- d[, .(data = sum(data, na.rm = TRUE)), by = .(name, year, unit, season, area)][data > 0]
    d[, dtype := lab][, keep, with = FALSE]
  }
  if(!is.null(lencomp)) out$len <- cmp(lencomp, "Length comps")
  if(!is.null(agecomp)) out$age <- cmp(agecomp, "Age comps")
  if(!length(out)) stop("no data given", call. = FALSE)
  d <- rbindlist(lapply(out, function(z) z[, name := as.character(name)]), use.names = TRUE)
  d[, rel := data/max(data, na.rm = TRUE), by = .(dtype, name)]
  d[, dtype := factor(dtype, levels = c("Catch", "Index", "Length comps", "Age comps"))]
  fl <- unique(d$name)
  d[, name := factor(name, levels = rev(fl))]
  d <- .add_group(d, dims)
  cols <- stats::setNames(ss3cols(length(fl)), fl)
  p <- ggplot2::ggplot(d, ggplot2::aes(year, name)) +
    ggplot2::geom_point(ggplot2::aes(size = rel, fill = name), shape = 21, colour = "grey20", stroke = 0.2) +
    ggplot2::scale_size_area(max_size = 3.5, guide = "none") +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::facet_grid(if(length(dims)) dtype ~ group else dtype ~ ., scales = "free_y", space = "free_y") +
    ggplot2::labs(x = "Year", y = NULL) + scale_x_year() + theme_ss3() +
    ggplot2::theme(strip.text.y = ggplot2::element_text(angle = 0))
  p
}
