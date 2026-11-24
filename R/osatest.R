# Tests of OSA residual distributions and model scorecards

#' Residuals by panel
#'
#' Long table of residuals with the panels of \code{\link{ss3plot_resid}()} (fleet and the 'by' dimensions; cells
#' pooled over dropped dimensions with Stouffer's sum, still standard normal for OSA residuals)
#'
#' @param x FLQuants of residuals (\code{\link{ss3_resid}()})
#' @param by dimensions kept as separate panels
#' @return data.table with columns panel, name, q, year and data
#' @noRd
.resid_panels <- function(x, by = "unit"){
  dims <- .bydims(by)
  d <- .flq_long(x, 1, "resid")[is.finite(data)]
  setnames(d, "quant", "q"); d[, q := as.numeric(as.character(q))]
  d <- d[, .(data = sum(data)/sqrt(.N)), by = c("name", "q", "year", dims)]
  if(!"unit" %in% dims) d[, unit := "unique"]
  d <- .add_group(d, dims); .panel(d, dims)
}

#' Anderson-Darling statistic against N(0,1)
#'
#' @param r residuals
#' @return A^2 for a fully specified standard normal
#' @noRd
.ad_stat <- function(r){
  r <- sort(r); n <- length(r)
  p <- pmin(pmax(stats::pnorm(r), 1e-12), 1 - 1e-12)
  -n - mean((2*seq_len(n) - 1)*(log(p) + log(1 - rev(p))))
}

#' Tests of OSA residual distributions
#'
#' Tests whether the OSA residuals of each composition (fleet and unit, as the panels of
#' \code{\link{ss3plot_resid}()}) follow their theoretical N(0,1) distribution: the Anderson-Darling statistic against a
#' fixed N(0,1) with a Monte Carlo p-value, and descriptive checks (mean near 0, SD near 1, the root mean square
#' distance of the QQ plot from the 1:1 line, and the share outside +-1.96, about 5\% expected).
#'
#' The Monte Carlo p-value assumes independent residuals; temporal or cohort dependence can make the nominal
#' p-values too small, and testing many fleets and models raises the chance of false failures. Treat the results as
#' flags, not as definitive judgements of model quality.
#'
#' @param x FLQuants of OSA residuals (\code{\link{ss3_resid}()}) or an SS3 report (residuals computed with
#'   \code{...})
#' @param by dimensions tested separately, any of \code{"unit"} (default), \code{"season"}, \code{"area"}
#' @param nsim Monte Carlo simulations for the p-value
#' @param alpha significance level (FAIL when p < alpha)
#' @param seed random seed
#' @param ... arguments to \code{\link{ss3_resid}()} when \code{x} is an SS3 report
#' @return data.table: panel, N, mean, SD, QQ RMSE, AD statistic, p-value, result (PASS/FAIL), \% outside +-1.96
#' @seealso \code{\link{ss3plot_resid}(test = TRUE)}, \code{\link{ss3_scorecard}()}
#' @examples
#' data(ss3simple)
#' r <- ss3_resid(ss3simple)
#' ss3_osatest(r, nsim = 199)
#' ss3plot_resid(r, test = TRUE)
#' ss3plot_resid(r, type = "qq", test = TRUE)
#' @export
ss3_osatest <- function(x, by = "unit", nsim = 999, alpha = 0.05, seed = 1, ...){
  if(.is_ss3rep(x)) x <- ss3_resid(x, ...)
  d <- .resid_panels(x, by)
  old <- if(exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if(!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed)
  out <- d[, {
    r <- data; n <- length(r); a <- .ad_stat(r)
    sim <- vapply(seq_len(nsim), function(i) .ad_stat(stats::rnorm(n)), 0)
    p <- (1 + sum(sim >= a))/(nsim + 1)
    .(N = n, mean = round(mean(r), 3), SD = round(stats::sd(r), 3),
      `QQ RMSE` = round(sqrt(mean((sort(r) - stats::qnorm(stats::ppoints(n)))^2)), 3),
      AD = round(a, 3), p = round(p, 3), result = ifelse(p < alpha, "FAIL", "PASS"),
      `outside 1.96 (%)` = round(100*mean(abs(r) > 1.96), 1))
  }, by = panel]
  setattr(out, "alpha", alpha); setattr(out, "nsim", nsim)
  out[]
}

#' Pointwise 95% envelopes of a normal QQ plot
#'
#' @param n sample size
#' @param level coverage
#' @return data.table of theoretical quantiles and the envelope of the ordered standard normal sample
#' @noRd
.qq_envelope <- function(n, level = 0.95){
  i <- seq_len(n); a <- (1 - level)/2
  data.table(theo = stats::qnorm(stats::ppoints(n)), lo = stats::qnorm(stats::qbeta(a, i, n - i + 1)),
             hi = stats::qnorm(stats::qbeta(1 - a, i, n - i + 1)))
}

#' Model diagnostic scorecard
#'
#' Pass/fail comparison of several models (Carvalho et al. 2021 diagnostic criteria): runs tests of the index
#' residuals (p >= \code{runs.p}), the hindcast MASE of each index (< \code{mase}), the retrospective and forecast
#' Mohn's rho of SSB (within \code{rho}), and the Anderson-Darling test of the OSA residuals of the compositions
#' (\code{\link{ss3_osatest}()}). Criteria that cannot be computed for a model (e.g. no retrospective runs) are NA.
#'
#' @param models named list of SS3 reports
#' @param retro named list (same names as \code{models}) of retrospective runs (lists of reports or
#'   \code{SSsummarize()} output, reference run first); models without are skipped for the retro and MASE criteria
#' @param quant composition data for the OSA test: \code{"len"}, \code{"age"} or \code{NULL} for none
#' @param rho interval of Mohn's rho of SSB that passes; default \code{c(-0.15, 0.2)} for long-lived species
#'   (Hurtado-Ferro et al. 2015; \code{c(-0.22, 0.3)} for short-lived species)
#' @param runs.p,mase thresholds of the runs test p-value and of MASE
#' @param nsim Monte Carlo simulations of the OSA test
#' @param ... further arguments to \code{\link{ss3_resid}()} (e.g. \code{nsamp = "input"})
#' @return list: \code{table} (wide, model by criterion, 1 pass, 0 fail, NA unavailable), \code{values} (long table
#'   of the statistics and results), \code{osa} (the OSA tests) and \code{mohn}
#' @seealso \code{\link{ss3plot_scorecard}()}
#' @examples
#' \donttest{
#' data(ss3sma); data(retro.sma)
#' sc <- ss3_scorecard(list(sma = ss3sma), retro = list(sma = retro.sma), quant = "len", nsim = 199)
#' sc$table
#' ss3plot_scorecard(sc)
#' }
#' @export
ss3_scorecard <- function(models, retro = NULL, quant = "len", rho = c(-0.15, 0.2), runs.p = 0.05, mase = 1,
                          nsim = 999, ...){
  if(.is_ss3rep(models)) models <- list(model = models)
  nm <- names(models) %||% paste0("model", seq_along(models)); names(models) <- nm
  vals <- list(); osa <- list(); mo <- list()
  add <- function(m, crit, value, pass) vals[[length(vals) + 1]] <<- data.table(model = m, criterion = crit,
                                                                                value = value, pass = as.integer(pass))
  for(m in nm){
    r <- models[[m]]
    ru <- tryCatch(ss3_runs(r), error = function(e) NULL)
    if(!is.null(ru)) for(i in seq_len(nrow(ru))) add(m, paste("Runs index:", ru$index[i]), ru$runs.p[i], ru$runs.p[i] >= runs.p)
    if(!is.null(retro[[m]])){
      rt <- tryCatch(ss3_retro(retro[[m]]), error = function(e) NULL)
      if(!is.null(rt)){
        mr <- ss3_mohn(rt, "SSB"); mo[[m]] <- cbind(model = m, ss3_mohn(rt, c("SSB", "F")))
        cb <- mr[peel == "Combined"]
        add(m, "Retro SSB rho", cb$rho, cb$rho >= rho[1] & cb$rho <= rho[2])
        add(m, "Forecast SSB rho", cb$fcrho, cb$fcrho >= rho[1] & cb$fcrho <= rho[2])
      }
      hc <- tryCatch(ss3_hindcast(retro[[m]])$mase, error = function(e) NULL)
      if(!is.null(hc)) for(i in seq_len(nrow(hc))) add(m, paste("MASE index:", hc$index[i]), hc$MASE[i],
                                                        if(is.finite(hc$MASE[i])) hc$MASE[i] < mase else NA)
    }
    if(!is.null(quant)){
      ot <- tryCatch(ss3_osatest(r, nsim = nsim, quant = quant, ...), error = function(e) NULL)
      if(!is.null(ot)){
        osa[[m]] <- cbind(model = m, ot)
        for(i in seq_len(nrow(ot))) add(m, paste0("OSA ", quant, ": ", ot$panel[i]), ot$p[i], ot$result[i] == "PASS")
      }
    }
  }
  v <- rbindlist(vals)
  lev <- unique(v$criterion)
  ord <- c(grep("^Runs", lev, value = TRUE), grep("rho$", lev, value = TRUE), grep("^MASE", lev, value = TRUE),
           grep("^OSA", lev, value = TRUE))
  w <- dcast(v, model ~ factor(criterion, levels = ord), value.var = "pass")
  w <- w[match(nm, model)]
  list(table = w[], values = v[], osa = rbindlist(osa), mohn = rbindlist(mo))
}

#' Plot a model diagnostic scorecard
#'
#' Tile plot of \code{\link{ss3_scorecard}()}: green pass, red fail, grey unavailable, with the statistic in each cell.
#'
#' @param x output of \code{\link{ss3_scorecard}()}
#' @param values print the statistic (p-value, rho, MASE) in each cell
#' @param text show the title, subtitle and caption (\code{FALSE}: none; labels inside the plot and legends stay)
#' @return ggplot
#' @examples
#' \donttest{
#' data(ss3sma); data(retro.sma)
#' ss3plot_scorecard(ss3_scorecard(list(sma = ss3sma), retro = list(sma = retro.sma), nsim = 199))
#' }
#' @export
ss3plot_scorecard <- function(x, values = TRUE, text = TRUE){
  v <- copy(x$values)
  crit <- setdiff(names(x$table), "model")
  v <- CJ(model = unique(v$model), criterion = crit)[v, on = c("model", "criterion"), `:=`(value = i.value, pass = i.pass)]
  v[, criterion := factor(criterion, levels = rev(crit))][, model := factor(model, levels = x$table$model)]
  v[, res := factor(ifelse(is.na(pass), "unavailable", ifelse(pass == 1, "pass", "fail")), levels = c("pass", "fail", "unavailable"))]
  v[, lab := ifelse(is.finite(value), formatC(value, digits = 2, format = "f"), "")]
  p <- ggplot2::ggplot(v, ggplot2::aes(model, criterion, fill = res)) + ggplot2::geom_tile(colour = "white") +
    ggplot2::scale_fill_manual(values = c(pass = "#a6d96a", fail = "#f4a582", unavailable = "grey85"), drop = FALSE,
                               name = NULL) +
    ggplot2::scale_x_discrete(position = "top") +
    ggplot2::labs(x = NULL, y = NULL, title = "Model diagnostics",
                  caption = "pass: runs test p >= 0.05; Mohn's rho of SSB within the interval;\nMASE < 1; OSA Anderson-Darling p >= 0.05") +
    theme_ss3() + ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "bottom")
  if(values) p <- p + ggplot2::geom_text(ggplot2::aes(label = lab), size = 3)
  .drop_text(p, text)
}

utils::globalVariables(c("panel", "res", "lab", "i.value", "i.pass", "criterion", "theo", "lo", "hi", "result"))
