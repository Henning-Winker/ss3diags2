# structure.R - selectivity, stock by area and stock-recruitment (r4ss HTML tabs Sel, Timeseries, S-R)

.sel_factors <- c(len = "Lsel", ret = "Ret", keep = "Keep", dead = "Dead", mort = "Mort", age = "Asel2")

#' Selectivity and retention by fleet
#'
#' Selectivity at length (\code{"len"}), retention (\code{"ret"}), discard mortality (\code{"mort"}), the kept
#' (\code{"keep"}) or dead (\code{"dead"}) fractions at length, or the derived selectivity at age (\code{"age"},
#' \code{Asel2}) from \code{ss3rep$sizeselex} / \code{ss3rep$ageselex}, as \code{FLQuants} by fleet with
#' len or age, year (the years in which SS3 reports a change, i.e. time blocks) and unit (sex) dimensions.
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param what \code{"len"}, \code{"ret"}, \code{"keep"}, \code{"dead"}, \code{"mort"} or \code{"age"}
#' @param fleets fleet numbers or names (default: all fleets)
#' @param years years to keep (default: the reported years from the start to the end year; blocks with identical
#'   selectivity are dropped)
#' @param long return the long data.table
#' @return FLQuants (named by fleet) or data.table
#' @export
ss3_sel <- function(ss3rep, what = c("len", "ret", "keep", "dead", "mort", "age"), fleets = NULL, years = NULL,
                    long = FALSE){
  what <- match.arg(what)
  quant <- if(what == "age") "age" else "len"
  s <- as.data.table(if(what == "age") ss3rep$ageselex else ss3rep$sizeselex)
  s <- s[Factor == .sel_factors[[what]]]
  if(!nrow(s)) stop("no ", .sel_factors[[what]], " in the report", call. = FALSE)
  if(what == "age") s <- s[Seas == min(Seas) & Morph %in% s[, min(Morph), by = Sex]$V1]
  bins <- names(s)[suppressWarnings(!is.na(as.numeric(names(s))))]
  fleets <- .pick_fleets(ss3rep, fleets, .fleetnames(ss3rep)[sort(unique(s$Fleet))])
  s <- s[Fleet %in% match(fleets, .fleetnames(ss3rep))]
  yrs <- years %||% sort(unique(s$Yr[s$Yr >= ss3rep$startyr & s$Yr <= ss3rep$endyr]))
  if(is.null(years)) yrs <- unique(c(ss3rep$startyr, yrs))
  s[Yr < ss3rep$startyr, Yr := ss3rep$startyr]
  s <- s[Yr %in% yrs][, .SD[.N], by = .(Fleet, Yr, Sex)]
  d <- melt(s, id.vars = c("Fleet", "Yr", "Sex"), measure.vars = bins, variable.name = "bin", value.name = "data",
            variable.factor = FALSE)
  # drop years identical to the previous reported year (no change)
  setorder(d, Fleet, Sex, bin, Yr)
  key <- d[, .(h = paste(signif(data, 6), collapse = ",")), by = .(Fleet, Yr)]
  key[, keep := c(TRUE, h[-1] != h[-.N]), by = Fleet]
  d <- merge(d, key[keep == TRUE, .(Fleet, Yr)], by = c("Fleet", "Yr"))
  nsx <- ss3rep$nsexes %||% 1
  d <- data.table(fleet = d$Fleet, name = .fleetnames(ss3rep)[d$Fleet], q = d$bin, year = d$Yr,
                  unit = .unit_lab(d$Sex, rep(1, nrow(d)), nsx), season = "all", area = "unique", iter = "1",
                  data = d$data)
  setnames(d, "q", quant)
  if(long) return(d[])
  long2flqs(d, quant = quant, units = "NA", order = fleets)
}

#' Plot selectivity by fleet
#'
#' Selectivity curves by fleet from \code{\link{ss3_sel}()}, one line per time block (colour) and unit (line type),
#' with an optional overlay of retention (dashed).
#'
#' @param sel FLQuants from \code{\link{ss3_sel}()}
#' @param ret optional FLQuants of retention (\code{ss3_sel(ss3rep, "ret")}), drawn dashed
#' @param ncol number of panel columns
#' @return ggplot
#' @export
ss3plot_sel <- function(sel, ret = NULL, ncol = NULL){
  prep <- function(x, lab){
    d <- .flq_long(x, 1, "sel")
    d[, x := as.numeric(as.character(quant))][, what := lab]
    d
  }
  d <- prep(sel, "selectivity")
  if(!is.null(ret)) d <- rbind(d, prep(ret[intersect(names(ret), names(sel))], "retention"))
  d[, name := factor(as.character(name), levels = names(sel))]
  d[, block := factor(year)]
  d[, lt := if(length(unique(unit)) > 1) paste(what, unit) else what]
  cols <- stats::setNames(ss3cols(nlevels(d$block)), levels(d$block))
  lv <- unique(d$lt[order(d$what != "selectivity", d$unit)])
  lts <- stats::setNames(c("solid", "dashed", "dotted", "dotdash")[seq_along(lv)], lv)
  ggplot2::ggplot(d, ggplot2::aes(x, data, colour = block, linetype = lt,
                                  group = interaction(block, lt))) +
    ggplot2::geom_line(linewidth = 0.7) + ggplot2::scale_colour_manual(values = cols) +
    ggplot2::scale_linetype_manual(values = lts, breaks = lv) +
    ggplot2::facet_wrap(~ name, ncol = ncol) +
    ggplot2::labs(x = if(is(sel[[1]], "FLQuant") && names(dimnames(sel[[1]]))[1] == "age") "Age" else "Length",
                  y = "Selectivity", colour = "From year", linetype = NULL) +
    theme_ss3() + ggplot2::theme(legend.title = ggplot2::element_text())
}

#' Spawning stock, biomass and recruitment by area
#'
#' Extracts the SS3 time series by area (\code{ss3rep$timeseries}): spawning output at spawning time, total biomass
#' at the beginning of the year and recruitment summed over seasons, as \code{FLQuants} with an area dimension.
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param years years to keep (default: model years)
#' @return FLQuants \code{SSB}, \code{B} and \code{Rec} (quant "all", area dimension)
#' @export
ss3_byarea <- function(ss3rep, years = NULL){
  ts <- as.data.table(ss3rep$timeseries)
  yrs <- years %||% seq(ss3rep$startyr, ss3rep$endyr)
  ts <- ts[Yr %in% yrs]
  na <- ss3rep$nareas %||% 1; sp <- ss3rep$spawnseas %||% 1
  mk <- function(x, val, u) long2flq(x[, .(age = "all", year = Yr, unit = "unique", season = "all",
                                           area = .area_lab(Area, na), iter = "1", data = get(val))],
                                     units = u, years = yrs)
  FLCore::FLQuants(SSB = mk(ts[Seas == sp], "SpawnBio", ss3rep$SpawnOutputUnits %||% "NA"),
                   B = mk(ts[Seas == min(Seas)], "Bio_all", "t"),
                   Rec = mk(ts[, .(Recruit_0 = sum(Recruit_0)), by = .(Yr, Area)], "Recruit_0", "1000"))
}

#' Plot stock quantities by area
#'
#' SSB, biomass and recruitment by area from \code{\link{ss3_byarea}()}, as stacked areas (\code{"area"}) or the
#' share of each area (\code{"share"}).
#'
#' @param x FLQuants from \code{\link{ss3_byarea}()}
#' @param type \code{"area"} (default) or \code{"share"}
#' @return ggplot
#' @export
ss3plot_byarea <- function(x, type = c("area", "share")){
  type <- match.arg(type)
  d <- .flq_long(x, 1, "x")
  d[, area := factor(paste("area", area), levels = paste("area", .lev(area)))]
  cols <- stats::setNames(ss3cols(nlevels(d$area)), levels(d$area))
  ggplot2::ggplot(d, ggplot2::aes(year, data, fill = area)) +
    ggplot2::geom_area(position = if(type == "share") "fill" else "stack", colour = "white", linewidth = 0.1) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::facet_wrap(~ name, scales = "free_y", ncol = 1, strip.position = "left") +
    ggplot2::labs(x = "Year", y = if(type == "share") "Share by area" else NULL, fill = NULL) +
    scale_x_year() + theme_ss3() + ggplot2::theme(strip.placement = "outside")
}

#' Plot the stock-recruitment relationship
#'
#' SS3 recruitment against spawning output (\code{ss3rep$recruit}) with the expected stock-recruit curve, coloured by
#' year, and the recruitment deviations over time.
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param type \code{"sr"} (stock-recruit, default) or \code{"devs"} (recruitment deviations)
#' @return ggplot
#' @export
ss3plot_sr <- function(ss3rep, type = c("sr", "devs")){
  type <- match.arg(type)
  rc <- as.data.table(ss3rep$recruit)[Yr >= ss3rep$startyr & Yr <= ss3rep$endyr]
  if(type == "devs")
    return(ggplot2::ggplot(rc[!is.na(dev)], ggplot2::aes(Yr, dev)) + ggplot2::geom_hline(yintercept = 0, colour = "grey50") +
             ggplot2::geom_segment(ggplot2::aes(xend = Yr, yend = 0), colour = "grey40") +
             ggplot2::geom_point(ggplot2::aes(colour = dev > 0), size = 2) +
             ggplot2::scale_colour_manual(values = c(`TRUE` = "#0000CC", `FALSE` = "#CC0000"), guide = "none") +
             ggplot2::labs(x = "Year", y = "Log recruitment deviation") + scale_x_year() + theme_ss3())
  sh <- .sr_shape(ss3rep, "ss3")
  p <- as.data.table(ss3rep$parameters)
  R0 <- exp(p[Label == "SR_LN(R0)"]$Value[1]); S0 <- ss3rep$SBzero
  s <- seq(0, max(c(S0, rc$SpawnBio), na.rm = TRUE), length.out = 200)
  cv <- data.table(SpawnBio = s, R = R0*sh(s)/sh(S0))
  ggplot2::ggplot(rc, ggplot2::aes(SpawnBio, pred_recr)) +
    ggplot2::geom_line(data = cv, ggplot2::aes(SpawnBio, R), colour = "black", linewidth = 0.8) +
    ggplot2::geom_path(colour = "grey70") +
    ggplot2::geom_point(ggplot2::aes(fill = Yr), shape = 21, size = 2.5) +
    ggplot2::scale_fill_gradientn(colours = ss3cols(8)) +
    ggplot2::expand_limits(x = 0, y = 0) +
    ggplot2::labs(x = paste0("Spawning output (", ss3rep$SpawnOutputUnits %||% "", ")"), y = "Recruitment",
                  fill = "Year") + theme_ss3() + ggplot2::theme(legend.title = ggplot2::element_text())
}
