# retro.R - retrospective analysis: FLStocks by peel, Mohn's rho, plot

#' Retrospective runs as FLStocks
#'
#' Converts a list of SS3 retrospective runs (\code{r4ss::SSgetoutput()} of \code{r4ss::retro()} runs, reference run
#' first) to \code{FLStocks} of \code{\link{ss3stk}()} MLE stocks, one per peel, named by the peel end year (taken
#' from \code{Retro_year}). With \code{fwd = TRUE} (default) each stock includes the year after its end year, for the
#' one-year-ahead (forecast) Mohn's rho and the forecast segments of \code{\link{ss3plot_retro}()}. Each stock carries
#' its end year in \code{attr(x, "endyr")}.
#'
#' @param retro list of SS3 reports (reference run first)
#' @param fwd include the year after each peel's end year (default TRUE)
#' @param ... passed to \code{ss3stk()}
#' @return FLStocks (named by peel end year)
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' rt <- ss3_retro(retro.sma)
#' ss3_mohn(rt)
#' ss3plot_retro(rt)
#' }
#' @export
ss3_retro <- function(retro, fwd = TRUE, ...){
  if(.is_ss3rep(retro)) stop("retro must be a list of SS3 reports (reference run first)", call. = FALSE)
  ey <- vapply(retro, .endyr, 0)
  stks <- lapply(seq_along(retro), function(i)
    suppressMessages(ss3stk(retro[[i]], output = "mle", fwd = fwd, name = as.character(ey[i]), ...)))
  out <- FLCore::FLStocks(stats::setNames(stks, as.character(ey)))
  out
}

# value of quantity q of a stock (median over iterations) by year
.qval <- function(stk, q){
  x <- ss3_flqs(stk, q)[[1]]
  stats::setNames(apply(x@.Data, 2, stats::median, na.rm = TRUE), dimnames(x)$year)
}

#' Mohn's rho
#'
#' Mohn's rho of retrospective peels against the reference run, \eqn{\rho_i = (x_i(T_i) - x_{ref}(T_i))/x_{ref}(T_i)}
#' at each peel's end year \eqn{T_i}, and the forecast (one-year-ahead) rho at \eqn{T_i + 1}, with their means
#' ("Combined"), as \code{ss3diags::SShcbias()}. Works for any FLStocks with the peels' end years (FLStocks from
#' \code{\link{ss3_retro}()}, or e.g. JABBA3 retrospective stocks), using \code{\link{ss3_flqs}()} quantities.
#'
#' @param x FLStocks, reference first (or a list of SS3 reports, passed to \code{ss3_retro()})
#' @param quants any quantities of \code{\link{ss3_flqs}()} (default "SSB", "F")
#' @param endyrs peel end years (default \code{attr(stk, "endyr")}, else the names of \code{x})
#' @return data.table: quant, peel, rho, fcrho (forecast rho; NA when the peel has no year after its end year)
#' @export
ss3_mohn <- function(x, quants = c("SSB", "F"), endyrs = NULL){
  if(!is(x, "FLStocks")) x <- ss3_retro(x)
  ey <- endyrs %||% vapply(seq_along(x), function(i) attr(x[[i]], "endyr") %||% as.numeric(names(x)[i]), 0)
  rbindlist(lapply(quants, function(q){
    ref <- .qval(x[[1]], q)
    r <- rbindlist(lapply(seq_along(x)[-1], function(i){
      v <- .qval(x[[i]], q); t <- as.character(ey[i]); t1 <- as.character(ey[i] + 1)
      data.table(quant = q, peel = t, rho = (v[t] - ref[t])/ref[t],
                 fcrho = if(t1 %in% names(v) && t1 %in% names(ref)) (v[t1] - ref[t1])/ref[t1] else NA_real_)
    }))
    rbind(r, data.table(quant = q, peel = "Combined", rho = mean(r$rho), fcrho = mean(r$fcrho)))
  }))
}

#' Plot retrospective analysis
#'
#' Retrospective peels of SSB, F, recruitment or stock status against the reference run with its confidence interval
#' (log-scale SEs of the reference stock), each peel to its end year and, with \code{fcast = TRUE}, the one-year-ahead
#' forecast as a dashed segment ending in a point; Mohn's rho (forecast rho in brackets) in the panel titles, as
#' \code{ss3diags::SSplotRetro()} and \code{JABBA3::jbplot_retro()}.
#'
#' @param x FLStocks from \code{\link{ss3_retro}()} (reference first), or a list of SS3 reports
#' @param type quantities of \code{\link{ss3_flqs}()} (default "SSB", "F")
#' @param fcast show the one-year-ahead forecasts
#' @param ci confidence level of the reference interval (0 for none)
#' @param xmin first year shown (default: 2 x number of peels before the last end year... all years if NULL)
#' @param ncol facet columns
#' @return ggplot
#' @export
ss3plot_retro <- function(x, type = c("SSB", "F"), fcast = TRUE, ci = 0.95, xmin = NULL, ncol = NULL){
  if(!is(x, "FLStocks")) x <- ss3_retro(x)
  ey <- vapply(seq_along(x), function(i) attr(x[[i]], "endyr") %||% as.numeric(names(x)[i]), 0)
  rho <- ss3_mohn(x, type, ey)[peel == "Combined"]
  lab <- stats::setNames(paste0(.ylab_quant(type), ": rho = ", sprintf("%.2f", rho$rho),
                                ifelse(fcast & is.finite(rho$fcrho), sprintf(" (%.2f)", rho$fcrho), "")), type)
  pl <- names(x)
  d <- rbindlist(lapply(seq_along(x), function(i) rbindlist(lapply(type, function(q){
    v <- .qval(x[[i]], q)
    data.table(peel = pl[i], quant = q, year = as.numeric(names(v)), est = v)
  }))))
  d[, endyr := ey[match(peel, pl)]]
  d[, peel := factor(peel, levels = pl)][, quant := factor(quant, levels = type)]
  ref <- d[peel == pl[1] & year <= ey[1]]
  # reference interval from the log-scale SEs (SSB, F, B/Btgt, F/Ftgt)
  se <- attr(x[[1]], "se.log"); z <- if(ci > 0) stats::qnorm(1 - (1 - ci)/2) else 0
  sel <- c(SSB = "SSB", F = "F", Rec = "Rec", `B/Btgt` = "stock", `F/Ftgt` = "harvest")
  ref[, c("lo", "hi") := NA_real_]
  if(!is.null(se) && z > 0) for(k in intersect(names(sel), type)){
    s <- as.numeric(se[[sel[[k]]]]); yy <- as.numeric(dimnames(se[[1]])$year)
    ref[quant == k, `:=`(lo = est*exp(-z*s[match(year, yy)]), hi = est*exp(z*s[match(year, yy)]))]
  }
  peels <- d[peel != pl[1] & year <= endyr]
  fc <- d[peel != pl[1] & year %in% c(endyr, endyr + 1)][order(peel, quant, year)]
  fc <- fc[, if(.N == 2) .SD, by = .(peel, quant)]
  xmin <- xmin %||% (min(ey) - max(10, 2*length(x)))
  cols <- stats::setNames(c("black", ss3cols(length(x) - 1)), pl)
  p <- ggplot2::ggplot(mapping = ggplot2::aes(year, est, colour = peel))
  if(any(is.finite(ref$lo))) p <- p + ggplot2::geom_ribbon(data = ref[is.finite(lo)], ggplot2::aes(ymin = lo, ymax = hi),
                                                           fill = "grey70", colour = NA, alpha = 0.5)
  p <- p + ggplot2::geom_line(data = ref, linewidth = 0.8) + ggplot2::geom_line(data = peels, linewidth = 0.6)
  if(fcast && nrow(fc)) p <- p + ggplot2::geom_line(data = fc, linetype = 2, linewidth = 0.5) +
    ggplot2::geom_point(data = fc[year == endyr + 1], ggplot2::aes(fill = peel), shape = 21, colour = "black", size = 1.8,
                        show.legend = FALSE) + ggplot2::scale_fill_manual(values = cols, guide = "none")
  p + ggplot2::scale_colour_manual(values = cols, breaks = pl, labels = c("Ref", pl[-1])) +
    ggplot2::facet_wrap(~ quant, scales = "free_y", ncol = ncol, labeller = ggplot2::as_labeller(lab)) +
    ggplot2::coord_cartesian(xlim = c(xmin, max(ey) + if(fcast) 1 else 0)) + ggplot2::expand_limits(y = 0) +
    ggplot2::labs(x = "Year", y = NULL) + scale_x_year() + theme_ss3()
}

utils::globalVariables(c("peel", "endyr", "rho", "fcrho"))
