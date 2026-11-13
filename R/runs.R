# runs.R - residual runs tests with 3-sigma limits, joint residuals

#' Index and mean length / age residuals
#'
#' Log residuals of the abundance indices, \eqn{\log(obs) - \log(exp)} (observations used in the SS3 likelihood), or of
#' the mean length / age by fleet and year (sexes and partitions pooled, as \code{ss3diags::SScompsTA1.8()}), as
#' \code{FLQuants} by fleet, the input of \code{\link{ss3_runs}()}, \code{\link{ss3plot_runs}()} and
#' \code{\link{ss3plot_jres}()}.
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param quant \code{"cpue"} (default), \code{"len"}, \code{"age"} (mean length / age), \code{"con"} (conditional
#'   age-at-length: mean age by length bin combined over bins weighted by sample size, as
#'   \code{ss3diags::SScompsTA1.8(type = "con")}) or \code{"size"} (generalised size compositions)
#' @return FLQuants (quant "all", year x unit x season x area)
#' @export
ss3_iresid <- function(ss3rep, quant = c("cpue", "len", "age", "con", "size")){
  quant <- match.arg(quant)
  d <- .hc_obs(ss3rep, quant)
  if(quant == "cpue") d <- d[use > 0 & !is.na(like)]
  ns <- ss3rep$nseasons %||% 1
  d <- d[, .(name, age = "all", year, unit = "unique", season = .season_lab(season, ns), area = "unique", iter = "1",
             data = log(obs) - log(exp))]
  d <- d[, .(data = sum(data)/sqrt(.N)), by = .(name, age, year, unit, season, area, iter)]
  out <- long2flqs(d, quant = "age", units = "log residual", order = unique(d$name))
  attr(out, "type") <- quant
  out
}

# Wald-Wolfowitz runs test on the signs (as randtests::runs.test(threshold = 0)); returns the p-value
.runs_p <- function(x, mixing = c("less", "two.sided", "greater")){
  mixing <- match.arg(mixing)
  s <- sign(x[is.finite(x)]); s <- s[s != 0]
  n1 <- sum(s > 0); n2 <- sum(s < 0); n <- n1 + n2
  if(n1 == 0 || n2 == 0) return(0.001)
  R <- 1 + sum(s[-1] != s[-n])
  mu <- 2*n1*n2/n + 1; v <- 2*n1*n2*(2*n1*n2 - n)/(n^2*(n - 1))
  z <- (R - mu)/sqrt(v)
  switch(mixing, less = stats::pnorm(z), greater = 1 - stats::pnorm(z), two.sided = 2*min(stats::pnorm(z), 1 - stats::pnorm(z)))
}

# 3-sigma limits from the average moving range (as ss3diags::ssruns_sig3)
.sig3 <- function(x, mu = 0){
  mr <- abs(diff(x - mu)); amr <- mean(mr, na.rm = TRUE)
  mr <- mr[mr < 3.267*amr]; amr <- mean(mr, na.rm = TRUE)
  stdev <- amr/1.128
  c(mu - 3*stdev, mu + 3*stdev)
}

#' Runs tests with three-sigma limits
#'
#' Runs test for randomness of the residual signs (Wald-Wolfowitz, one-sided for positive autocorrelation by default)
#' and three-sigma control limits from the average moving range, by index, as \code{ss3diags::SSrunstest()}
#' (Carvalho et al. 2021, Fish. Res. 240: 105959; Winker et al. 2018). Series with 3 or fewer residuals or spanning 3
#' years or less are excluded.
#'
#' @param x FLQuants of residuals (e.g. \code{\link{ss3_iresid}()}, OSA residuals) or an SS3 report (index residuals)
#' @param mixing \code{"less"} (default: too few runs, positive autocorrelation), \code{"two.sided"} or \code{"greater"}
#' @param ... passed to \code{ss3_iresid()} when \code{x} is an SS3 report (e.g. \code{quant = "len"})
#' @return data.table: index, runs.p, test (Passed, Failed, Excluded), sigma3.lo, sigma3.hi
#' @export
ss3_runs <- function(x, mixing = "less", ...){
  if(.is_ss3rep(x)) x <- ss3_iresid(x, ...)
  d <- .flq_long(x, 1, "resid")[is.finite(data)][order(name, year)]
  rbindlist(lapply(levels(d$name), function(f){
    r <- d[name == f]$data; yy <- d[name == f]$year
    if(length(r) <= 3 || diff(range(yy)) <= 3)
      return(data.table(index = f, runs.p = NA_real_, test = "Excluded", sigma3.lo = NA_real_, sigma3.hi = NA_real_))
    p <- round(.runs_p(r, mixing), 3); lim <- .sig3(r)
    data.table(index = f, runs.p = p, test = ifelse(p < 0.05, "Failed", "Passed"), sigma3.lo = lim[1], sigma3.hi = lim[2])
  }))
}

#' Plot runs tests
#'
#' Residuals by index with the three-sigma limits shaded green (runs test passed) or red (failed), residuals beyond the
#' limits in red, and the runs test p-value in the panel title, as \code{ss3diags::SSplotRunstest()} and
#' \code{ggplotFL::plotRunstest()}.
#'
#' @param x FLQuants of residuals or an SS3 report
#' @param mixing runs test alternative, see \code{\link{ss3_runs}()}
#' @param ncol facet columns
#' @param ... passed to \code{ss3_iresid()} when \code{x} is an SS3 report
#' @return ggplot
#' @export
ss3plot_runs <- function(x, mixing = "less", ncol = NULL, ...){
  if(.is_ss3rep(x)) x <- ss3_iresid(x, ...)
  rt <- ss3_runs(x, mixing)
  d <- .flq_long(x, 1, "resid")[is.finite(data)]
  d[, name := as.character(name)]
  d <- merge(d, rt, by.x = "name", by.y = "index")[test != "Excluded"]
  lab <- rt[test != "Excluded", .(index, lab = paste0(index, ": p = ", sprintf("%.3f", runs.p)))]
  d[, panel := factor(lab$lab[match(name, lab$index)], levels = lab$lab)]
  d[, out := data < sigma3.lo | data > sigma3.hi]
  bx <- d[, .(xmin = min(year) - 0.5, xmax = max(year) + 0.5, lo = sigma3.lo[1], hi = sigma3.hi[1], test = test[1]),
          by = panel]
  lim <- d[, .(y = c(min(-1, sigma3.lo[1]*1.3, min(data)), max(1, sigma3.hi[1]*1.3, max(data))), year = year[1]), by = panel]
  ylab <- switch(attr(x, "type") %||% "cpue", cpue = "Index residuals", len = "Mean length residuals",
                 age = "Mean age residuals", con = "Mean age-at-length residuals", size = "Mean size residuals", "Residuals")
  ggplot2::ggplot(d, ggplot2::aes(year, data)) +
    ggplot2::geom_rect(data = bx, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = lo, ymax = hi, fill = test),
                       inherit.aes = FALSE, alpha = 0.6) +
    ggplot2::scale_fill_manual(values = c(Passed = "#cbe368", Failed = "#ef8575"), name = "Runs test") +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_segment(ggplot2::aes(xend = year, y = 0, yend = data)) +
    ggplot2::geom_point(ggplot2::aes(colour = out), size = 2) +
    ggplot2::geom_point(colour = "black", shape = 21, size = 2) +
    ggplot2::scale_colour_manual(values = c(`FALSE` = "white", `TRUE` = "#d50102"), guide = "none") +
    ggplot2::geom_blank(data = lim, ggplot2::aes(year, y)) +
    ggplot2::facet_wrap(~ panel, scales = "free", ncol = ncol) +
    ggplot2::labs(x = "Year", y = ylab) + scale_x_year(n = 4) + theme_ss3()
}

#' Plot joint residuals
#'
#' Residuals of all indices (or mean lengths) together by year, with boxplots across indices, a loess smoother of
#' the mean residual and the overall RMSE (in percent) in the title, as \code{ss3diags::SSplotJABBAres()}.
#'
#' @param x FLQuants of residuals or an SS3 report
#' @param ... passed to \code{ss3_iresid()} when \code{x} is an SS3 report
#' @return ggplot; the RMSE table (by index and combined) in \code{attr(p, "rmse")}
#' @export
ss3plot_jres <- function(x, ...){
  if(.is_ss3rep(x)) x <- ss3_iresid(x, ...)
  d <- .flq_long(x, 1, "resid")[is.finite(data)]
  d[, data := ifelse(abs(data) > 3, NA_real_, data)]
  rmse <- rbind(d[, .(index = as.character(name[1]), RMSE.perc = round(100*sqrt(mean(data^2, na.rm = TRUE)), 1),
                      nobs = sum(is.finite(data))), by = name][, name := NULL],
                data.table(index = "Combined", RMSE.perc = round(100*sqrt(mean(d$data^2, na.rm = TRUE)), 1),
                           nobs = sum(is.finite(d$data))))
  nl <- nlevels(d$name)
  off <- if(nl > 1) stats::setNames(seq(-0.25, 0.25, length.out = nl), levels(d$name)) else stats::setNames(0, levels(d$name))
  d[, x := year + off[as.character(name)]]
  mn <- d[, .(m = mean(data, na.rm = TRUE)), by = year][order(year)]
  mn[, sm := if(.N > 4) stats::predict(stats::loess(m ~ year)) else m]
  cols <- stats::setNames(ss3cols(nl), levels(d$name))
  ylab <- switch(attr(x, "type") %||% "cpue", cpue = "Index residuals", len = "Mean length residuals",
                 age = "Mean age residuals", con = "Mean age-at-length residuals", size = "Mean size residuals", "Residuals")
  p <- ggplot2::ggplot(d) +
    ggplot2::geom_boxplot(ggplot2::aes(year, data, group = year), fill = "grey85", colour = "grey50", outlier.shape = NA,
                          width = 0.7, na.rm = TRUE) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) +
    ggplot2::geom_segment(ggplot2::aes(x = x, xend = x, y = 0, yend = data, colour = name), na.rm = TRUE) +
    ggplot2::geom_point(ggplot2::aes(x, data, fill = name), shape = 21, size = 2, na.rm = TRUE) +
    ggplot2::geom_line(data = mn, ggplot2::aes(year, sm), linewidth = 1) +
    ggplot2::scale_colour_manual(values = cols, guide = "none") + ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(x = "Year", y = ylab, title = paste0("RMSE = ", rmse[index == "Combined"]$RMSE.perc, "%")) +
    scale_x_year() + theme_ss3()
  attr(p, "rmse") <- rmse
  p
}

utils::globalVariables(c("test", "runs.p", "sigma3.lo", "sigma3.hi", "out", "m", "sm", "RMSE.perc"))
