# ss3plot.R - one-page summary of an SS3 assessment (as JABBA3's jbplot summary)

#' Observed indices scaled to SSB
#'
#' Observed indices scaled to SSB (geometric mean ratio over the index years), for the SSB and B/Btgt panels
#'
#' @param ss3rep SS3 report (\code{r4ss::SS_output()})
#' @param ssb SSB by year
#' @return data.table of indices scaled by the geometric mean ratio over the index years
#' @noRd
.idx2ssb <- function(ss3rep, ssb){
  idx <- tryCatch(ss3_indices(ss3rep, long = TRUE), error = function(e) NULL)
  if(is.null(idx) || !nrow(idx)) return(NULL)
  idx <- idx[is.finite(obs) & obs > 0, .(obs = mean(obs)), by = .(name, year)]
  idx <- merge(idx, ssb, by = "year")
  idx[, sc := exp(mean(log(ssb) - log(obs))), by = name]
  idx[, .(name, year, value = obs*sc)]
}

#' Recruitment deviations with 95% intervals
#'
#' Recruitment deviations with 95% intervals from the SS3 parameter SEs
#'
#' @param ss3rep SS3 report (\code{r4ss::SS_output()})
#' @return data.table year, dev, lower, upper
#' @noRd
.recdevs <- function(ss3rep){
  rp <- as.data.table(ss3rep$recruitpars)
  if(!nrow(rp)) return(NULL)
  rp <- rp[Yr >= ss3rep$startyr & Yr <= .endyr(ss3rep)]
  rp[, .(year = Yr, dev = Value, lo = Value - 1.96*Parm_StDev, hi = Value + 1.96*Parm_StDev,
         type = sub("_RecrDev", "", type))]
}

#' Panel strip label for a ggplot
#'
#' @param p ggplot
#' @param lab label
#' @return ggplot
#' @noRd
.strip <- function(p, lab) p + ggplot2::facet_wrap(ggplot2::vars(!!lab)) +
  ggplot2::theme(plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank())

#' One-page summary plot of an SS3 assessment
#'
#' A 3 x 3 summary of an SS3 assessment in the style of the JABBA3 summary plot: SSB, F and catch with reference
#' points (top), B/Btgt and F/Ftgt with the Kobe plot (middle), and the stock-recruitment curve, the recruitment
#' deviations and the production phase plot (bottom). Biomass is SSB. Uncertainty comes from MVLN draws
#' (\code{\link{ss3stk}()}); the observed indices are shown in the SSB and B/Btgt panels scaled to SSB (geometric mean
#' ratio over the years of each index), as a visual check of the trend only.
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param stk optional FLStockR with iterations from \code{ss3stk(ss3rep, output = "iters")} (built if \code{NULL})
#' @param iters number of MVLN draws when \code{stk} is built
#' @param lag Kobe lag: 0 pairs F(y) with B(y), 1 pairs F(y) with B(y + 1) (see \code{\link{ss3plot_kobe}()})
#' @param indices show the observed indices scaled to SSB
#' @param legend position of the index legend, drawn once outside the panels: \code{"bottom"} (default),
#'   \code{"right"}, \code{"top"} or \code{"none"}; \code{TRUE} / \code{FALSE} for \code{"bottom"} / \code{"none"}
#' @param fwd show the SS3 forecast years in the trajectories (default \code{FALSE}: assessment years only)
#' @param probs quantiles for the uncertainty bands
#' @param text show the title, subtitle and caption (\code{FALSE}: none; labels inside the plot and legends stay)
#' @return a patchwork of ggplots (needs the patchwork package)
#' @seealso \code{\link{ss3plot_trj}()}, \code{\link{ss3plot_kobe}()}, \code{\link{ss3plot_sr}()},
#'   \code{\link{ss3plot_production}()}
#' @examples
#' pra <- readRDS(system.file("extdata", "pra.27.3a4a", "pra.27.3a4a.rds", package = "ss3diags2"))
#' ss3plot(pra, iters = 200, lag = 1)
#' @export
ss3plot <- function(ss3rep, stk = NULL, iters = 1000, lag = 0, indices = TRUE, fwd = FALSE,
                    probs = c(0.05, 0.2, 0.5, 0.8, 0.95), legend = c("bottom", "right", "top", "none"), text = TRUE){
  .drop_text((function(){
  if(is.logical(legend)) legend <- if(isTRUE(legend[1])) "bottom" else "none"
  legend <- match.arg(legend, c("bottom", "right", "top", "none"))
  .need("patchwork", "package 'patchwork' is required for ss3plot()")
  if(!.is_ss3rep(ss3rep)) stop("ss3rep must be an r4ss::SS_output() list", call. = FALSE)
  stk <- stk %||% ss3stk(ss3rep, output = "iters", iters = iters, fwd = TRUE)
  ey <- attr(stk, "endyr") %||% .endyr(ss3rep)
  st <- if(fwd) stk else FLCore::window(stk, end = ey)
  attributes(st)[c("basis", "endyr", "mle", "se.log", "eqcurve")] <- attributes(stk)[c("basis", "endyr", "mle", "se.log", "eqcurve")]
  noleg <- ggplot2::theme(legend.position = "none")
  trj <- function(q) ss3plot_trj(st, type = q, probs = probs) + ggplot2::labs(x = NULL, caption = NULL) + noleg
  pB <- trj("SSB"); pBr <- trj("B/Btgt")
  if(indices){
    ssb <- data.table(year = as.numeric(dimnames(stk)$year), ssb = c(FLCore::iterMedians(FLCore::ssb(stk))))
    ix <- .idx2ssb(ss3rep, ssb)
    if(!is.null(ix) && nrow(ix)){
      bt <- .rp(stk, "Btgt")
      cols <- stats::setNames(ss3cols(length(unique(ix$name))), unique(ix$name))
      pt <- function(d) list(ggplot2::geom_point(data = d, ggplot2::aes(year, value, fill = name), shape = 21,
                                                 size = 1.8, inherit.aes = FALSE),
                             ggplot2::scale_fill_manual(values = cols, name = NULL),
                             ggplot2::guides(fill = ggplot2::guide_legend(nrow = if(legend %in% c("bottom", "top"))
                               ceiling(length(cols)/5) else NULL, override.aes = list(size = 2.5))))
      # the index legend is drawn once, collected outside the panels (see legend)
      pB <- pB + pt(ix) + ggplot2::theme(legend.position = legend)
      if(is.finite(bt)) pBr <- pBr + pt(ix[, .(name, year, value = value/bt)]) + noleg
    }
  }
  pK <- .strip(ss3plot_kobe(stk, lag = lag) + noleg, "Kobe")
  pSR <- .strip(ss3plot_sr(ss3rep) + noleg + ggplot2::labs(x = "SSB"), "Stock-recruitment")
  rd <- .recdevs(ss3rep)
  pR <- if(is.null(rd)) ggplot2::ggplot() + ggplot2::theme_void() else
    .strip(ggplot2::ggplot(rd, ggplot2::aes(year, dev)) +
             ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = "grey45", alpha = 0.3) +
             ggplot2::geom_hline(yintercept = 0, linetype = 2) + ggplot2::geom_line(linewidth = 0.6) +
             ggplot2::labs(x = NULL, y = "log scale") + scale_x_year() + theme_ss3(), "Recruitment deviation")
  pP <- .strip(ss3plot_production(ss3stk(ss3rep), type = "phase") + noleg + ggplot2::labs(x = "SSB"),
               "Surplus production phase")
  patchwork::wrap_plots(pB, trj("F"), trj("Catch"), pBr, trj("F/Ftgt"), pK, pSR, pR, pP, ncol = 3) +
    patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(theme = ggplot2::theme(legend.position = legend))
  })(), text)
}
