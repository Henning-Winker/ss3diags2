# plots_resid.R - composition residual plots (FLR inputs), growth / cohort lines as FLicc

# one facet per fleet x group with data, free x (years with data only)
.panel <- function(d, dims){
  d[, panel := if(length(dims)) paste(name, group, sep = " | ") else as.character(name)]
  lv <- unique(d[order(name, group)]$panel)
  d[, panel := factor(panel, levels = lv)]
  d
}

# cohort tracks: von Bertalanffy length (or age) by year for birth years every 'by' years, per panel
.cohort_lines <- function(d, quant, lh, by = 2, timing = 0.5){
  pan <- d[, .(unit = unit[1], ymin = min(year), ymax = max(year), qmin = min(q), qmax = max(q)), by = panel]
  rbindlist(lapply(seq_len(nrow(pan)), function(i){
    z <- pan[i]
    z$unit <- as.character(z$unit)
    z$unit <- as.character(z$unit)
    units <- if(z$unit == "unique" && !is.null(lh) && "unit" %in% names(dimnames(lh))) dimnames(lh)$unit else z$unit
    rbindlist(lapply(units, function(u){
      yy <- seq(z$ymin - 0.5, z$ymax + 0.5, by = 0.1)
      if(quant == "age"){
        births <- seq(floor(z$ymin - z$qmax), ceiling(z$ymax - z$qmin), by = by)
        out <- CJ(birth = births, year = yy)[, q := year - birth]
      } else {
        vb <- .vb_unit(lh, u)
        if(!all(is.finite(vb))) return(NULL)
        amin <- max(vb[["t0"]] - log(1 - min(z$qmin/vb[["linf"]], 0.99))/vb[["k"]], 0)
        amax <- vb[["t0"]] - log(1 - min(0.95, z$qmax/vb[["linf"]]))/vb[["k"]]
        births <- seq(floor(z$ymin - amax), ceiling(z$ymax - amin), by = by)
        out <- CJ(birth = births, year = yy)[, q := vb[["linf"]]*(1 - exp(-vb[["k"]]*(year - birth + timing - vb[["t0"]])))]
      }
      out[year >= birth & q >= z$qmin & q <= z$qmax][, `:=`(panel = z$panel, lunit = u)]
    }))
  }))
}

#' Plot composition residuals
#'
#' Residuals of the length or age compositions from \code{FLQuants} (e.g. \code{\link{ss3_resid}()}: OSA or Pearson),
#' as in FLicc: bubbles by year and length (or age), one panel per fleet and unit (e.g. sex) with data, only the years
#' with data, blue positive, red negative; cohort tracks (von Bertalanffy growth of year classes from \code{lhpar},
#' e.g. \code{\link{ss3_lhpar}()}, by unit; diagonals for ages) to recognise year classes growing through the
#' compositions. \code{type = "len"} pools the residuals over years by bin (Stouffer, \eqn{\sum r/\sqrt{k}}, still
#' standard normal for OSA residuals) to show persistent misfit at given lengths (selectivity, Linf); \code{"qq"}
#' gives normal QQ plots by fleet. OSA residuals are conditional on the smaller bins: too few small fish gives negative
#' residuals at small lengths and positive ones at larger lengths (the conditional share of the larger bins is higher),
#' so read patterns along length as conditional, not as a direct excess or deficit of large fish.
#'
#' @param x FLQuants of residuals (len or age x year x unit x season x area)
#' @param type \code{"bubble"} (default), \code{"len"} (pooled over years by bin) or \code{"qq"}
#' @param by dimensions kept as separate panels, any of \code{"unit"} (default), \code{"season"}, \code{"area"}
#' @param lhpar FLPar with \code{linf}, \code{k} (and \code{t0}), with or without a unit dimension, for the cohort
#'   tracks of length compositions (not needed for ages); \code{NULL} for none
#' @param cohorts draw cohort tracks (default TRUE when possible)
#' @param cohort_by years between cohort tracks
#' @param cap residuals beyond +-cap drawn at the cap (bubble and len plots)
#' @param ncol number of panel columns
#' @return ggplot
#' @examples
#' data(ss3ank)
#' r <- ss3_resid(ss3ank)
#' ss3plot_resid(r, lhpar = ss3_lhpar(ss3ank))
#' ss3plot_resid(r, type = "len")
#' ss3plot_resid(r, type = "qq")
#' @export
ss3plot_resid <- function(x, type = c("bubble", "len", "qq"), by = "unit", lhpar = NULL, cohorts = TRUE,
                          cohort_by = 2, cap = 4, ncol = NULL){
  type <- match.arg(type); dims <- .bydims(by)
  rtype <- attr(x, "type") %||% "osa"
  lab <- paste(if(rtype == "osa") "OSA" else "Pearson", "residual")
  qn <- FLCore::quant(if(is(x, "FLQuant")) x else x[[1]])
  d <- .flq_long(x, 1, "resid")[is.finite(data)]
  setnames(d, "quant", "q"); d[, q := as.numeric(as.character(q))]
  # several cells pooled when 'by' drops dimensions: Stouffer
  keep <- c("name", "q", "year", dims)
  d <- d[, .(data = sum(data)/sqrt(.N)), by = keep]
  if(!"unit" %in% dims) d[, unit := "unique"]
  d <- .add_group(d, dims); d <- .panel(d, dims)
  xlab <- switch(qn, len = "Length", age = "Age", qn)
  if(type == "qq"){
    return(ggplot2::ggplot(d, ggplot2::aes(sample = data)) +
             ggplot2::stat_qq(size = 0.7, alpha = 0.5) + ggplot2::stat_qq_line(colour = "#D6604D") +
             ggplot2::facet_wrap(~ panel, ncol = ncol) +
             ggplot2::labs(x = "Theoretical N(0,1) quantiles", y = paste(lab, "quantiles")) + theme_ss3())
  }
  d[, r := pmax(pmin(data, cap), -cap)]
  d[, sign := factor(ifelse(data >= 0, "positive", "negative"), levels = c("positive", "negative"))]
  fills <- c(positive = "#2166AC", negative = "#D6604D")
  if(type == "len"){
    a <- d[, .(data = sum(data)/sqrt(.N)), by = .(panel, q)]
    a[, `:=`(r = pmax(pmin(data, cap), -cap), capped = abs(data) > cap)]
    return(ggplot2::ggplot(a, ggplot2::aes(q, r)) +
             ggplot2::geom_hline(yintercept = c(-2, 2), linetype = "dashed", colour = "grey60") +
             ggplot2::geom_hline(yintercept = 0, colour = "grey30") +
             ggplot2::geom_col(ggplot2::aes(fill = ifelse(r >= 0, "positive", "negative")), width = .bw(a$q)) +
             ggplot2::scale_fill_manual(values = fills, guide = "none") +
             ggplot2::facet_wrap(~ panel, ncol = ncol, scales = "free_x") +
             ggplot2::labs(x = xlab, y = paste(lab, "(years pooled)")) + theme_ss3())
  }
  p <- ggplot2::ggplot(d, ggplot2::aes(year, q))
  if(isTRUE(cohorts) && (qn == "age" || !is.null(lhpar))){
    tr <- .cohort_lines(d, qn, lhpar, cohort_by)
    if(nrow(tr)) p <- p + ggplot2::geom_line(data = tr, ggplot2::aes(year, q, group = interaction(birth, lunit),
                                                                     linetype = lunit),
                                             colour = "grey55", linewidth = 0.3) +
      ggplot2::scale_linetype_manual(values = c(F = 1, M = 2, unique = 1), guide = if(length(unique(tr$lunit)) > 1)
        ggplot2::guide_legend(title = "growth") else "none")
  }
  p + ggplot2::geom_point(ggplot2::aes(size = abs(r), fill = sign), shape = 21, colour = "grey25", stroke = 0.15,
                          alpha = 0.85) +
    ggplot2::scale_fill_manual(values = fills, labels = c(positive = "obs > exp", negative = "obs < exp")) +
    ggplot2::scale_size_area(max_size = 3.5, limits = c(0, cap), name = paste0("|", lab, "|")) +
    ggplot2::facet_wrap(~ panel, scales = "free", ncol = ncol) +
    ggplot2::labs(x = "Year", y = xlab) + scale_x_year(n = 4) + theme_ss3() +
    ggplot2::theme(legend.title = ggplot2::element_text(size = 8))
}

utils::globalVariables(c("panel", "birth", "lunit", "q", "r", "sign", "capped"))
