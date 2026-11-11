# plots_aux.R - theme, r4ss colours and helpers shared by the ss3plot_ functions

#' ggplot2 theme for ss3diags2
#'
#' The JABBA3 theme (\code{theme_bw()} with light strips, small legend keys, no minor grid), so plots from
#' ss3diags2, JABBA3 and FLicc can be combined.
#' @param base_size base font size
#' @param legend.position legend position (default "bottom")
#' @export
theme_ss3 <- function(base_size = 11, legend.position = "bottom"){
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(legend.key.size = ggplot2::unit(0.3, "cm"), legend.key.height = ggplot2::unit(0.5, "cm"),
                   legend.key.width = ggplot2::unit(0.5, "cm"), legend.title = ggplot2::element_blank(),
                   legend.text = ggplot2::element_text(size = 8), legend.position = legend.position,
                   strip.background = ggplot2::element_rect(fill = "grey92"),
                   panel.spacing.x = ggplot2::unit(0.8, "lines"), panel.grid.minor = ggplot2::element_blank())
}

#' r4ss colours
#'
#' Colour conventions of \code{r4ss}, so ss3diags2 plots look like the familiar SS3 plots:
#' \itemize{
#'   \item \code{ss3cols(n)}: fleets and other categories, \code{r4ss::rich.colors.short()} (first colour dropped
#'     for n > 2, as \code{r4ss::SSplotCatch()})
#'   \item \code{ss3seascols(n)}: seasons, black for one season, blue4, green3, orange2, red3 for four
#'     (\code{r4ss::SSplotIndices()})
#'   \item \code{ss3unitcols}: units (sexes) as r4ss composition fits: combined ("unique") green, F red, M blue
#'   \item \code{ss3bubblecols}: bubbles, positive red, negative blue, observations grey (\code{r4ss::SSplotComps()})
#' }
#' @param n number of colours
#' @param alpha transparency
#' @examples
#' ss3cols(5); ss3seascols(4); ss3unitcols
#' @export
ss3cols <- function(n, alpha = 1){
  rich <- function(n){
    x <- seq(0, 1, length.out = n); r <- 1/(1 + exp(20 - 35*x)); g <- pmin(pmax(0, -0.8 + 6*x - 5*x^2), 1)
    b <- stats::dnorm(x, 0.25, 0.15)/max(stats::dnorm(x, 0.25, 0.15))
    grDevices::rgb(r, g, b, alpha = alpha)
  }
  if(n > 2) rich(n + 1)[-1] else rich(n)
}

#' @rdname ss3cols
#' @export
ss3seascols <- function(n, alpha = 1){
  col <- if(n == 1) "black" else if(n == 4) c("blue4", "green3", "orange2", "red3") else ss3cols(n)
  if(alpha < 1) col <- grDevices::adjustcolor(col, alpha.f = alpha)
  col
}

#' @rdname ss3cols
#' @export
ss3unitcols <- c(unique = "#008000", F = "#CC0000", M = "#0000CC")

#' @rdname ss3cols
#' @export
ss3bubblecols <- c(positive = "#FF0000B3", negative = "#0000FFB3", obs = "#1A1A1AB3")

#' Integer year axis
#' @param n approximate number of breaks
#' @param ... passed to \code{scale_x_continuous()}
#' @export
scale_x_year <- function(n = 5, ...) ggplot2::scale_x_continuous(breaks = function(x) unique(round(pretty(x, n = n))), ...)

# --- FLR inputs -----------------------------------------------------------------------------------------------

# FLQuant / FLQuants / list of FLQuant -> long data.table with column 'name' (one iteration)
.flq_long <- function(x, iter = 1, what = "x"){
  if(is(x, "FLQuant")) x <- FLCore::FLQuants(stats::setNames(list(x), what))
  if(is.list(x) && !is(x, "FLQuants")) x <- FLCore::FLQuants(x)
  if(!is(x, "FLQuants")) stop("expected FLQuant or FLQuants", call. = FALSE)
  if(is.null(names(x)) || any(names(x) == "")) names(x) <- paste0(what, seq_along(x))
  x <- FLCore::FLQuants(lapply(x, function(q) if(dim(q)[6] > 1) FLCore::iter(q, iter) else q))
  d <- flqs2long(x)
  setnames(d, names(d)[1], "quant")
  d[, name := factor(name, levels = names(x))]
  d[]
}

# FLIndices / FLIndex(Biomass) -> long table: index (obs), index.var (se), timing
.fli_long <- function(x, iter = 1){
  if(is(x, "FLIndex") || is(x, "FLIndexBiomass")) x <- FLCore::FLIndices(stats::setNames(list(x), FLCore::name(x)))
  if(!is(x, "FLIndices")) stop("expected FLIndices or FLIndexBiomass", call. = FALSE)
  nm <- names(x); if(is.null(nm) || any(nm == "")) nm <- vapply(x, FLCore::name, "")
  rbindlist(lapply(seq_along(x), function(k){
    i <- x[[k]]
    q <- FLCore::index(i); v <- FLCore::index.var(i)
    if(dim(q)[6] > 1) q <- FLCore::iter(q, iter)
    if(dim(v)[6] > 1) v <- FLCore::iter(v, iter)
    d <- as.data.table(methods::as(q, "data.frame")); setnames(d, names(d)[1], "quant")
    dv <- as.data.table(methods::as(v, "data.frame"))
    d[, se := dv$data[match(paste(d$quant, year, unit, season, area), paste(dv[[1]], dv$year, dv$unit, dv$season,
                                                                            dv$area))]]
    tf <- i@range[["startf"]]; if(!is.finite(tf)) tf <- 0.5
    d[, `:=`(name = nm[k], timef = tf, units = FLCore::units(q))]
    d
  }), use.names = TRUE)[!is.na(data)][, name := factor(name, levels = nm)][]
}

# --- grouping by unit / season / area -------------------------------------------------------------------------

.bydims <- function(by){
  if(is.null(by) || !length(by)) return(character(0))
  bad <- setdiff(by, c("unit", "season", "area"))
  if(length(bad)) stop("'by' must be any of \"unit\", \"season\", \"area\"; got ", paste(bad, collapse = ", "),
                       call. = FALSE)
  unique(by)
}

# sum over the dimensions not in 'keep' (numbers, catch); 'value' columns are summed
.collapse <- function(d, keep, value = "data"){
  drop <- setdiff(c("unit", "season", "area"), keep)
  by <- setdiff(names(d), c(value, drop, "iter"))
  d <- d[, lapply(.SD, function(v) if(all(is.na(v))) NA_real_ else sum(v, na.rm = TRUE)), by = by, .SDcols = value]
  for(k in drop) d[, (k) := switch(k, unit = "unique", season = "all", area = "unique")]
  d
}

# 'group' column combining the 'by' dimensions, with readable labels and an ordered factor
.add_group <- function(d, dims){
  if(!length(dims)){ d[, group := factor("all")]; return(d) }
  lab1 <- function(k, v) switch(k, unit = as.character(v),
                                season = ifelse(v == "all", "all", paste0("season ", v)),
                                area = ifelse(v == "unique", "unique", paste0("area ", v)))
  lab <- do.call(paste, c(lapply(dims, function(k) lab1(k, d[[k]])), sep = " | "))
  ord <- do.call(order, lapply(dims, function(k) match(as.character(d[[k]]), .lev(d[[k]]))))
  d[, group := factor(lab, levels = unique(lab[ord]))]
  d
}

# r4ss colours for the groups
.group_cols <- function(d, dims){
  lv <- levels(d$group)
  if(!length(dims)) return(stats::setNames("black", lv))
  if(identical(dims, "unit")) {
    u <- as.character(d$unit[match(lv, d$group)])
    return(stats::setNames(ifelse(u %in% names(ss3unitcols), ss3unitcols[u], ss3cols(length(lv))), lv))
  }
  if(identical(dims, "season")) return(stats::setNames(ss3seascols(length(lv)), lv))
  stats::setNames(ss3cols(length(lv)), lv)
}

utils::globalVariables(c("group", "prop", "pexp", "mu", "mexp", "lo", "hi", "size", "dtype", "rel", "value",
                         "variable", "fit", "quant", "timef", "t", "xval", "ntot", "bw", "yfac"))
