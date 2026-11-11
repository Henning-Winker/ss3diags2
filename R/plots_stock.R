# plots_stock.R - trajectories, Kobe and production function from FLStockR (as JABBA3 jbplot_trj, jbplot_kobe,
# jbplot_spphase / jbplot_spdyn)

.kobe_fills <- c(green = "#2ca25f", yellow = "#fee391", orange = "#fd8d3c", red = "#de2d26")

# FLStockR / FLStock / FLStocks / list -> named list of stocks
.stklist <- function(x){
  if(is(x, "FLStock")) return(stats::setNames(list(x), FLCore::name(x) %||% "stock"))
  if(is(x, "FLStocks") || is.list(x)){
    nm <- names(x); if(is.null(nm) || any(nm == "")) nm <- vapply(x, FLCore::name, "")
    if(any(duplicated(nm)) || any(nm == "")) nm <- paste0("run", seq_along(x))
    return(stats::setNames(as.list(x), nm))
  }
  stop("expected FLStockR, FLStock, FLStocks or a list of them", call. = FALSE)
}

.rp <- function(stk, p){
  if(!is(stk, "FLStockR")) return(NA_real_)
  r <- FLCore::refpts(stk)
  if(p %in% dimnames(r)$params) stats::median(as.numeric(r[p, ])) else NA_real_
}

# long table of quantities for one stock: year, quant, iter, data
.stk_quants <- function(stk, type){
  B <- FLCore::ssb(stk); F <- FLCore::fbar(stk); R <- FLCore::rec(stk); C <- FLCore::catch(stk)
  btgt <- .rp(stk, "Btgt"); ftgt <- .rp(stk, "Ftgt"); msy <- .rp(stk, "MSY"); b0 <- .rp(stk, "B0")
  st <- attr(stk, "status")
  q <- list(SSB = B, F = F, Rec = R, Catch = C, `B/Btgt` = if(!is.null(st)) st$stock else B/btgt,
            `F/Ftgt` = if(!is.null(st)) st$harvest else F/ftgt, `B/B0` = B/b0, `Catch/MSY` = C/msy)
  bad <- setdiff(type, names(q)); if(length(bad)) stop("unknown type: ", paste(bad, collapse = ", "), call. = FALSE)
  rbindlist(lapply(type, function(k){
    d <- as.data.table(methods::as(q[[k]], "data.frame"))
    d[, .(year, iter, data, quant = k)]
  }))
}

# summary: median and quantiles from iterations, else MLE with lognormal intervals from attr(stk, "se.log")
.stk_summary <- function(stk, type, probs){
  d <- .stk_quants(stk, type)
  ni <- length(unique(d$iter))
  if(ni > 1){
    s <- d[, as.list(stats::quantile(data, probs, na.rm = TRUE, names = FALSE)), by = .(year, quant)]
    setnames(s, paste0("V", seq_along(probs)), c("lci", "lo1", "est", "hi1", "uci"))
    mle <- attr(stk, "mle")
    if(!is.null(mle)){                                   # MVLN draws: line at the MLE, as SSdeltaMVLN
      b0 <- .rp(stk, "B0"); msy <- .rp(stk, "MSY")
      m <- rbindlist(list(mle[, .(year, quant = "SSB", data = SSB)], mle[, .(year, quant = "F", data = F)],
                          mle[, .(year, quant = "Rec", data = Recr)], mle[, .(year, quant = "Catch", data = Catch)],
                          mle[, .(year, quant = "B/Btgt", data = stock)], mle[, .(year, quant = "F/Ftgt", data = harvest)],
                          mle[, .(year, quant = "B/B0", data = SSB/b0)], mle[, .(year, quant = "Catch/MSY", data = Catch/msy)]))
      s[m, est := i.data, on = .(year, quant)]
    }
    return(s)
  }
  se <- attr(stk, "se.log")
  z <- stats::qnorm(probs)
  sel <- c(SSB = "SSB", F = "F", Rec = "Rec", `B/Btgt` = "stock", `F/Ftgt` = "harvest", `B/B0` = "SSB")
  d[, c("lci", "lo1", "hi1", "uci") := NA_real_]
  if(!is.null(se)) for(k in intersect(names(sel), type)){
    s <- as.numeric(se[[sel[[k]]]]); yrs <- as.numeric(dimnames(se[[1]])$year)
    sv <- s[match(d[quant == k]$year, yrs)]
    d[quant == k, `:=`(lci = data*exp(z[1]*sv), lo1 = data*exp(z[2]*sv), hi1 = data*exp(z[4]*sv),
                       uci = data*exp(z[5]*sv))]
  }
  setnames(d, "data", "est"); d[, iter := NULL]; d
}

.ylab_quant <- function(type){
  c(SSB = "SSB", F = "F", Rec = "Recruitment", Catch = "Catch", `B/Btgt` = "B/Btgt", `F/Ftgt` = "F/Ftgt",
    `B/B0` = "B/B0", `Catch/MSY` = "Catch/MSY")[type]
}

#' Plot stock trajectories
#'
#' Trajectories of SSB, F, recruitment, catch and stock status (B/Btgt, F/Ftgt, B/B0, Catch/MSY) with nested
#' 90% and 60% bands (\code{probs}), as \code{JABBA3::jbplot_trj()} and \code{FLRef::plotAdvice()}, from an
#' \code{FLStockR} (e.g. \code{\link{ss3_FLStockR}()}, \code{FLRef::ss2FLStockR()}, \code{JABBA3::jb3_FLStockR()},
#' \code{FLicc::flicc2FLStockR()}). Bands come from the iterations (MVLN or MCMC), or for a single iteration from the
#' log-scale SEs in \code{attr(x, "se.log")}. Several stocks (a list or \code{FLStocks}) are coloured by run, so an
#' SS3 run compares directly with JABBA3 or FLicc. Reference points from \code{refpts()} are drawn as labelled lines.
#'
#' @param x FLStockR, FLStocks or a named list of FLStockR
#' @param type quantities, any of "SSB", "F", "Rec", "Catch", "B/Btgt", "F/Ftgt", "B/B0", "Catch/MSY"
#' @param probs five probabilities: outer band, inner band, median, inner, outer
#' @param refpts draw reference points (several runs: relative quantities only)
#' @param colours run colours (several runs; default \code{ss3cols()})
#' @param colour band colour for a single run
#' @param ncol facet columns
#' @param label.size size of the reference point labels
#' @return ggplot
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' ss3plot_trj(ss3_FLStockR(ss3sma))
#' ss3plot_trj(ss3_FLStockR(ss3sma, iters = 1000), type = c("B/Btgt", "F/Ftgt", "Catch"), ncol = 1)
#' }
#' @export
ss3plot_trj <- function(x, type = c("SSB", "F", "Rec", "Catch", "B/Btgt", "F/Ftgt"),
                        probs = c(0.05, 0.2, 0.5, 0.8, 0.95), refpts = TRUE, colours = NULL, colour = "grey45",
                        ncol = NULL, label.size = 2.5){
  sl <- .stklist(x); nr <- length(sl)
  d <- rbindlist(lapply(names(sl), function(n) .stk_summary(sl[[n]], type, probs)[, run := n]), fill = TRUE)
  miss <- d[, .(ok = any(is.finite(est))), by = quant][ok == FALSE]$quant
  if(length(miss)){
    message("no values (reference points missing?) for: ", paste(miss, collapse = ", "))
    type <- setdiff(type, miss); d <- d[quant %in% type]
  }
  d[, quant := factor(quant, levels = type)][, run := factor(run, levels = names(sl))]
  dr <- d[is.finite(lci) & is.finite(uci)]
  p <- ggplot2::ggplot(d, ggplot2::aes(year, est))
  if(nr == 1){
    if(nrow(dr)) p <- p + ggplot2::geom_ribbon(data = dr, ggplot2::aes(ymin = lci, ymax = uci), fill = colour, alpha = 0.25) +
      ggplot2::geom_ribbon(data = dr, ggplot2::aes(ymin = lo1, ymax = hi1), fill = colour, alpha = 0.45)
    p <- p + ggplot2::geom_line(linewidth = 0.6)
  } else {
    cols <- colours %||% ss3cols(nr)
    cols <- if(!is.null(names(cols)) && all(names(sl) %in% names(cols))) cols[names(sl)] else
      stats::setNames(rep_len(cols, nr), names(sl))
    if(nrow(dr)) p <- p + ggplot2::geom_ribbon(data = dr, ggplot2::aes(ymin = lci, ymax = uci, fill = run), alpha = 0.15) +
      ggplot2::geom_ribbon(data = dr, ggplot2::aes(ymin = lo1, ymax = hi1, fill = run), alpha = 0.25)
    p <- p + ggplot2::geom_line(ggplot2::aes(colour = run), linewidth = 0.7) +
      ggplot2::scale_colour_manual(values = cols) + ggplot2::scale_fill_manual(values = cols)
  }
  p <- p + ggplot2::facet_wrap(~ quant, scales = "free_y", ncol = ncol,
                               labeller = ggplot2::as_labeller(.ylab_quant(type))) +
    ggplot2::expand_limits(y = 0) + ggplot2::labs(x = "Year", y = NULL) + scale_x_year() + theme_ss3()
  if(refpts){
    s1 <- sl[[1]]
    rel <- data.table(quant = c("B/Btgt", "F/Ftgt", "B/B0", "Catch/MSY", "B/Btgt", "B/Btgt"),
                      name = c("Btgt", "Ftgt", "Btgt", "MSY", "Blim", "Bpa"),
                      value = c(1, 1, .rp(s1, "Btgt")/.rp(s1, "B0"), 1, .rp(s1, "Blim")/.rp(s1, "Btgt"),
                                .rp(s1, "Bpa")/.rp(s1, "Btgt")))
    abs <- data.table(quant = c("SSB", "SSB", "SSB", "F", "Catch", "Rec"), name = c("Btgt", "Blim", "Bpa", "Ftgt", "MSY", "R0"),
                      value = c(.rp(s1, "Btgt"), .rp(s1, "Blim"), .rp(s1, "Bpa"), .rp(s1, "Ftgt"), .rp(s1, "MSY"), .rp(s1, "R0")))
    refs <- if(nr == 1) rbind(rel, abs) else rel[name != "Btgt" | quant != "B/B0"]
    refs <- refs[quant %in% type & is.finite(value)]
    if(nrow(refs)){
      refs[, quant := factor(quant, levels = type)]
      lty <- c(Btgt = 2, Ftgt = 2, MSY = 2, R0 = 3, Blim = 3, Bpa = 4)
      p <- p + ggplot2::geom_hline(data = refs, ggplot2::aes(yintercept = value, linetype = name), colour = "darkred",
                                   linewidth = 0.5, show.legend = FALSE) +
        ggplot2::scale_linetype_manual(values = lty) +
        ggplot2::geom_text(data = refs, ggplot2::aes(x = min(d$year), y = value, label = name), hjust = 0, vjust = -0.4,
                           size = label.size, colour = "darkred", inherit.aes = FALSE)
    }
  }
  if(nr == 1) p <- p + ggplot2::theme(legend.position = "none")
  p
}

#' Kobe plot
#'
#' B/Btgt against F/Ftgt: the trajectory (MLE, or median of the iterations) and the joint uncertainty in one year from
#' the iterations (MVLN or MCMC) with highest-density contours and the percentage of draws in each quadrant, as
#' \code{JABBA3::jbplot_kobe()} and \code{ss3diags::SSplotKobe()}.
#'
#' @param x FLStockR with iterations (e.g. \code{ss3_FLStockR(ss3rep, iters = 1000)}); with one iteration only the
#'   trajectory is drawn
#' @param year year of the uncertainty cloud (default last year)
#' @param npoints maximum number of points shown
#' @param contours highest-density contour levels (\code{NULL} or \code{FALSE} for none; needs MASS)
#' @param style \code{"points"} (draws with contour lines) or \code{"contours"} (filled 50/80/95% regions)
#' @param labels label years on the trajectory every \code{labels} years (0 for first and last only)
#' @return ggplot
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' ss3plot_kobe(ss3_FLStockR(ss3sma, iters = 1000))
#' }
#' @export
ss3plot_kobe <- function(x, year = NULL, npoints = 2000, contours = c(0.5, 0.8, 0.95),
                         style = c("points", "contours"), labels = 0){
  style <- match.arg(style)
  if(!is(x, "FLStockR")) stop("x must be an FLStockR with refpts Btgt and Ftgt", call. = FALSE)
  btgt <- .rp(x, "Btgt"); ftgt <- .rp(x, "Ftgt")
  if(!is.finite(btgt) || !is.finite(ftgt)) stop("Kobe plot needs refpts Btgt and Ftgt", call. = FALSE)
  st <- attr(x, "status")
  bb <- if(!is.null(st)) st$stock else FLCore::ssb(x)/btgt
  ff <- if(!is.null(st)) st$harvest else FLCore::fbar(x)/ftgt
  yrs <- as.numeric(dimnames(bb)$year)
  mle <- attr(x, "mle")
  trj <- if(!is.null(mle)) data.table(year = mle$year, stock = mle$stock, harvest = mle$harvest) else
    data.table(year = yrs, stock = apply(bb@.Data, 2, stats::median, na.rm = TRUE),
               harvest = apply(ff@.Data, 2, stats::median, na.rm = TRUE))
  y <- year %||% max(trj$year)
  trj <- trj[year <= y & is.finite(stock) & is.finite(harvest)]
  pt <- if(dim(bb)[6] > 1){
    j <- match(as.character(y), dimnames(bb)$year)
    data.table(stock = as.numeric(bb[, j]), harvest = as.numeric(ff[, j]))
  } else NULL
  qd <- data.table(xmin = c(1, 0, 1, 0), xmax = c(Inf, 1, Inf, 1), ymin = c(0, 0, 1, 1), ymax = c(1, 1, Inf, Inf),
                   zone = c("green", "yellow", "orange", "red"))
  if(!is.null(pt)){
    zone <- with(pt, ifelse(stock >= 1 & harvest <= 1, "green", ifelse(stock < 1 & harvest <= 1, "yellow",
                 ifelse(stock >= 1 & harvest > 1, "orange", "red"))))
    pct <- round(100*prop.table(table(factor(zone, qd$zone))), 1)
    qd[, lab := sprintf("%.1f%%", pct[zone])]
    if(nrow(pt) > npoints) pt <- pt[sample(.N, npoints)]
  } else qd[, lab := zone]
  xmax <- max(2, trj$stock, if(!is.null(pt)) stats::quantile(pt$stock, 0.995))*1.05
  ymax <- max(2, trj$harvest, if(!is.null(pt)) stats::quantile(pt$harvest, 0.995), na.rm = TRUE)*1.05
  ends <- trj[c(1, .N)]
  lab <- if(labels > 0) trj[year %% labels == 0] else ends
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = qd, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = zone), alpha = 0.85) +
    ggplot2::scale_fill_manual(values = .kobe_fills, breaks = qd$zone, labels = qd$lab,
                               name = if(!is.null(pt)) paste("Uncertainty", y) else NULL,
                               guide = if(is.null(pt)) "none" else "legend")
  if(!is.null(pt)){
    if(style == "points") p <- p + ggplot2::geom_point(data = pt, ggplot2::aes(stock, harvest), colour = "grey15",
                                                       alpha = 0.3, size = 0.7)
    if(!isFALSE(contours) && length(contours)) p <- p + .kobe_hdr(pt, contours, filled = style == "contours",
                                                                 lims = c(0, xmax, 0, ymax))
  }
  p + ggplot2::geom_path(data = trj, ggplot2::aes(stock, harvest), linewidth = 0.6) +
    ggplot2::geom_point(data = trj, ggplot2::aes(stock, harvest), shape = 21, fill = "white", size = 1.3) +
    ggplot2::geom_point(data = ends, ggplot2::aes(stock, harvest), shape = c(22, 24), fill = "white", size = 3) +
    ggplot2::geom_text(data = lab, ggplot2::aes(stock, harvest, label = year), vjust = -1.1, size = 3) +
    ggplot2::coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax), expand = FALSE) +
    ggplot2::labs(x = paste0("B/", attr(x, "basis")$B %||% "Btgt"), y = paste0("F/", attr(x, "basis")$F %||% "Ftgt")) +
    theme_ss3(legend.position = "right") + ggplot2::theme(legend.title = ggplot2::element_text(size = 8))
}

# highest-density regions of the (stock, harvest) cloud: contour lines or filled 50/80/95% regions
.kobe_hdr <- function(pt, probs, filled = FALSE, lims){
  if(!requireNamespace("MASS", quietly = TRUE) || nrow(pt) < 50) return(NULL)
  k <- MASS::kde2d(pt$stock, pt$harvest, n = 120,
                   lims = c(range(pt$stock)*c(0.7, 1.3), range(pt$harvest)*c(0.7, 1.3)))
  z <- sort(k$z, decreasing = TRUE); cz <- cumsum(z)/sum(z)
  lev <- vapply(probs, function(p) z[which(cz >= p)[1]], 0)
  g <- data.table(expand.grid(stock = k$x, harvest = k$y)); g[, z := as.vector(k$z)]
  if(!filled) return(ggplot2::geom_contour(data = g, ggplot2::aes(stock, harvest, z = z), breaks = sort(lev),
                                           colour = "black", linewidth = 0.4))
  cfill <- c("cornsilk4", "grey", "cornsilk2")[seq_along(probs)]
  o <- order(probs, decreasing = TRUE)
  lapply(seq_along(o), function(i){
    L <- grDevices::contourLines(k$x, k$y, k$z, levels = lev[o[i]])
    if(!length(L)) return(NULL)
    cl <- rbindlist(lapply(seq_along(L), function(j) data.table(grp = j, x = L[[j]]$x, y = L[[j]]$y)))
    ggplot2::geom_polygon(data = cl, ggplot2::aes(x, y, group = grp), fill = cfill[i], colour = NA, alpha = 0.85)
  })
}

#' Production function plots
#'
#' The SS3 equilibrium yield curve against SSB with the catch trajectory (\code{type = "phase"}, as
#' \code{JABBA3::jbplot_spphase()}: zones of the Kobe plot in yield-biomass space), or the realised surplus production
#' SP(t) = SSB(t+1) - SSB(t) + C(t) against SSB with the equilibrium curve (\code{type = "dyn"}, as
#' \code{JABBA3::jbplot_spdyn()}).
#'
#' @param x FLStockR (MLE, or the median of iterations is used), e.g. \code{\link{ss3_FLStockR}()}
#' @param eq equilibrium curve: data.frame with columns \code{SSB} and \code{yield} (default
#'   \code{attr(x, "eqcurve")}, set by \code{ss3_FLStockR()} from \code{\link{ss3_eqcurve}()})
#' @param type \code{"phase"} (default) or \code{"dyn"}
#' @param rel relative axes: SSB/B0 and yield/MSY
#' @param spm overlay the Pella-Tomlinson fit to the equilibrium curve, extrapolated to SSB = 0 and K: \code{TRUE}
#'   (default: fitted with \code{\link{ss3_spm}()}), an \code{ss3spm} object, or \code{FALSE}
#' @return ggplot
#' @examples
#' \dontrun{
#' data(natl.sma, package = "ss3diags")
#' stk <- ss3_FLStockR(ss3sma)
#' ss3plot_production(stk)
#' ss3plot_production(stk, type = "dyn")
#' }
#' @export
ss3plot_production <- function(x, eq = NULL, type = c("phase", "dyn"), rel = FALSE, spm = TRUE){
  type <- match.arg(type)
  if(!is(x, "FLStockR")) stop("x must be an FLStockR", call. = FALSE)
  if(identical(attr(x, "basis")$ssbunits, "numbers"))
    warning("SSB is spawning output in numbers, not biomass: yield and surplus production are not on the same scale",
            call. = FALSE)
  eq <- as.data.table(eq %||% attr(x, "eqcurve"))
  if(!nrow(eq)) stop("no equilibrium curve: give 'eq' or build x with ss3_FLStockR()", call. = FALSE)
  med <- function(q) apply(q@.Data, 2, stats::median, na.rm = TRUE)
  B <- med(FLCore::ssb(x)); C <- med(FLCore::catch(x)); yrs <- as.numeric(dimnames(FLCore::ssb(x))$year)
  msy <- .rp(x, "MSY"); bmsy <- .rp(x, "Btgt"); b0 <- .rp(x, "B0")
  if(is.finite(attr(eq, "Bmsy") %||% NA)) bmsy <- attr(eq, "Bmsy")
  if(is.finite(attr(eq, "MSY") %||% NA)) msy <- attr(eq, "MSY")
  sx <- if(rel) b0 else 1; sy <- if(rel) msy else 1
  trj <- data.table(year = yrs, B = B/sx, C = C/sy, SP = (c(B[-1], NA) - B + C)/sy)
  cur <- data.table(b = eq$SSB/sx, sp = eq$yield/sy, what = "SS3 equilibrium")
  sp <- if(inherits(spm, "ss3spm")) spm else if(isTRUE(spm)) tryCatch(suppressWarnings(ss3_spm(eq)), error = function(e) NULL)
  if(!is.null(sp)) cur <- rbind(cur, data.table(b = sp$curve$SSB/sx, sp = sp$curve$yield/sy, what = "Pella-Tomlinson"))
  cur[, what := factor(what, levels = c("SS3 equilibrium", "Pella-Tomlinson"))]
  ltys <- c(`SS3 equilibrium` = 1, `Pella-Tomlinson` = 2); lws <- c(`SS3 equilibrium` = 1.1, `Pella-Tomlinson` = 0.6)
  Bm <- bmsy/sx; M <- msy/sy
  sel <- trj[unique(c(1, round(0.7*.N), .N))][, lab := factor(year, levels = year)]
  xl <- if(rel) "SSB/B0" else "SSB"
  if(type == "phase"){
    xmax <- max(cur$b, trj$B, na.rm = TRUE)*1.02; ymax <- max(trj$C*1.05, M*1.15, na.rm = TRUE); big <- 10*max(xmax, ymax)
    H <- M/Bm
    zones <- rbind(data.table(zone = "red", x = c(0, Bm, Bm, 0), y = c(0, M, big, big)),
                   data.table(zone = "yellow", x = c(0, Bm, Bm), y = c(0, 0, M)),
                   data.table(zone = "green", x = c(Bm, big, big, Bm), y = c(0, 0, H*big, M)),
                   data.table(zone = "orange", x = c(Bm, big, big, Bm), y = c(M, H*big, big, big)))
    p <- ggplot2::ggplot() +
      ggplot2::geom_polygon(data = zones, ggplot2::aes(x, y, group = zone, fill = zone), alpha = 0.85) +
      ggplot2::scale_fill_manual(values = .kobe_fills, guide = "none") +
      ggplot2::geom_line(data = cur, ggplot2::aes(b, sp, linetype = what, linewidth = what), colour = "blue") +
      ggplot2::geom_hline(yintercept = M, colour = "blue", linetype = 5) +
      ggplot2::annotate("segment", x = Bm, xend = Bm, y = 0, yend = M, colour = "blue", linetype = 2) +
      ggplot2::geom_path(data = trj, ggplot2::aes(B, C), na.rm = TRUE) +
      ggplot2::geom_point(data = trj, ggplot2::aes(B, C), size = 1.1, na.rm = TRUE) +
      ggplot2::geom_point(data = sel, ggplot2::aes(B, C, shape = lab), fill = "white", size = 3.2) +
      ggplot2::coord_cartesian(xlim = c(0, xmax), ylim = c(0, ymax), expand = FALSE) +
      ggplot2::labs(x = xl, y = if(rel) "Catch/MSY" else "Yield")
  } else {
    p <- ggplot2::ggplot() +
      ggplot2::geom_hline(yintercept = 0, linetype = 2) +
      ggplot2::geom_hline(yintercept = M, linetype = 3, colour = "grey40") +
      ggplot2::geom_vline(xintercept = Bm, linetype = 3, colour = "grey40") +
      ggplot2::geom_line(data = cur, ggplot2::aes(b, sp, linetype = what, linewidth = what), colour = "blue") +
      ggplot2::geom_path(data = trj, ggplot2::aes(B, SP), na.rm = TRUE) +
      ggplot2::geom_point(data = trj, ggplot2::aes(B, SP), size = 1.1, na.rm = TRUE) +
      ggplot2::geom_point(data = sel[is.finite(SP)], ggplot2::aes(B, SP, shape = lab), fill = "white", size = 3.2) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(c(0, 0.03)), limits = c(0, NA)) +
      ggplot2::labs(x = xl, y = if(rel) "SP/MSY" else "Surplus production",
                    subtitle = "realised SP(t) = SSB(t+1) - SSB(t) + C(t)")
  }
  p + ggplot2::scale_shape_manual(values = c(22, 21, 24), name = NULL) +
    ggplot2::scale_linetype_manual(values = ltys, name = NULL) + ggplot2::scale_linewidth_manual(values = lws, name = NULL) +
    theme_ss3()
}

utils::globalVariables(c("est", "lci", "uci", "lo1", "hi1", "run", "i.data", "xmin", "xmax", "ymin", "ymax", "zone",
                         "lab", "B", "C", "SP", "b", "sp", "grp", "z", "x", "y", "yield", "what", "ok"))
