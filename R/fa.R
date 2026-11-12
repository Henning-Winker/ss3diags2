# fa.R - fishing mortality at age by fleet, partial F, F by area and partial impact

.fa_id <- c("Area", "Fleet", "Sex", "Morph", "Yr", "Seas", "Era")

#' Fishing mortality at age by fleet
#'
#' Extracts F-at-age by fleet from \code{ss3rep$fatage} as \code{FLQuants}, one \code{FLQuant} per catch fleet with
#' the full SS3 structure: age, year, unit (sex), season and area (the area of the fleet). F-at-age is the
#' instantaneous rate on an annual scale within each season (apical F times selectivity), as in the SS3 report; the
#' annual F of a fleet is the sum over seasons weighted by the season durations (\code{ss3rep$seasdurations}).
#' Platoons / morphs within a sex are averaged.
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param fleets catch fleet numbers or names (default: all catch fleets)
#' @param era SS3 era(s) to keep: \code{"TIME"} (default), \code{"FORE"}, \code{"INIT"}
#' @param annual sum the seasons weighted by duration (season dimension "all")
#' @param long return the long data.table instead of FLQuants
#' @return FLQuants (named by fleet, units "f") or a data.table with columns \code{fleet}, \code{name}, \code{age},
#'   \code{year}, \code{unit}, \code{season}, \code{area}, \code{iter}, \code{data}
#' @references Code adapted from getFa() in the GSA 17-18 hake advice (H. Winker).
#' @examples
#' \dontrun{
#' fa <- ss3_fa(ss3rep)
#' ss3plot_fa(fa)
#' }
#' @export
ss3_fa <- function(ss3rep, fleets = NULL, era = "TIME", annual = FALSE, long = FALSE){
  fa <- as.data.table(ss3rep$fatage)
  if(!nrow(fa)) stop("ss3rep$fatage is empty (F_at_age not in the report)", call. = FALSE)
  fl <- ss3_fleets(ss3rep)
  fleets <- .pick_fleets(ss3rep, fleets, fl[type %in% c("catch", "bycatch") & fleet %in% fa$Fleet]$name)
  ages <- setdiff(names(fa), .fa_id)
  fa <- fa[Era %in% era & Fleet %in% match(fleets, .fleetnames(ss3rep))]
  d <- melt(fa, id.vars = intersect(.fa_id, names(fa)), measure.vars = ages, variable.name = "age",
            value.name = "data", variable.factor = FALSE)
  d <- d[, .(data = mean(data)), by = .(Fleet, Area, Sex, Yr, Seas, age)]
  ns <- ss3rep$nseasons %||% 1; na <- ss3rep$nareas %||% 1; nsx <- ss3rep$nsexes %||% 1
  if(annual){
    dur <- ss3rep$seasdurations %||% rep(1/ns, ns)
    d <- d[, .(data = sum(data*dur[Seas])), by = .(Fleet, Area, Sex, Yr, age)][, Seas := 1L]
    ns <- 1
  }
  d <- data.table(fleet = d$Fleet, name = .fleetnames(ss3rep)[d$Fleet], age = d$age, year = d$Yr,
                  unit = .unit_lab(d$Sex, rep(1, nrow(d)), nsx), season = .season_lab(d$Seas, ns),
                  area = .area_lab(d$Area, na), iter = "1", data = d$data)
  if(long) return(d[])
  long2flqs(d, quant = "age", units = "f", order = fleets)
}

# ages of the SS3 F_std basis ("..._for_ages_unweighted_1_3"); NULL if not age-based
.fbar_ages <- function(ss3rep){
  b <- ss3rep$F_std_basis %||% ss3rep$F_report_basis %||% ""
  m <- regmatches(b, regexec("ages_[a-z_]*?([0-9]+)_([0-9]+)", b))[[1]]
  if(length(m) == 3) as.numeric(m[2]):as.numeric(m[3]) else NULL
}

# fleet groups: NULL (each fleet), "area" (fleets by area) or a named list of fleet names / numbers
.fleet_groups <- function(ss3rep, groups, fleets){
  nm <- .fleetnames(ss3rep)
  if(is.null(groups)) return(stats::setNames(as.list(fleets), fleets))
  if(identical(groups, "area")){
    fl <- ss3_fleets(ss3rep)[name %in% fleets]
    return(split(fl$name, paste0("area ", fl$area)))
  }
  if(!is.list(groups) || is.null(names(groups))) stop("'groups' must be NULL, \"area\" or a named list of fleets",
                                                      call. = FALSE)
  groups <- lapply(groups, function(g) if(is.numeric(g)) nm[g] else g)
  bad <- setdiff(unlist(groups), fleets)
  if(length(bad)) stop("unknown catch fleet(s) in 'groups': ", paste(bad, collapse = ", "), call. = FALSE)
  miss <- setdiff(fleets, unlist(groups))
  if(length(miss)) warning("fleet(s) not in any group (ignored): ", paste(miss, collapse = ", "), call. = FALSE)
  groups
}

# catch at age in numbers by fleet, summed over seasons, sexes and areas: data.table fleet, year, age, data
.catage_annual <- function(ss3rep){
  ca <- as.data.table(ss3rep$catage)
  ca <- ca[, !duplicated(names(ca)), with = FALSE]
  if("Era" %in% names(ca)) ca <- ca[Era == "TIME"]
  if("Type" %in% names(ca)) ca <- ca[Type == "dead" | is.na(Type)]
  ages <- names(ca)[grepl("^[0-9]+$", names(ca))]
  d <- melt(ca, id.vars = c("Fleet", "Yr"), measure.vars = ages, variable.name = "age", value.name = "data",
            variable.factor = FALSE)
  d[, .(data = sum(data)), by = .(fleet = Fleet, year = Yr, age = as.numeric(age))]
}

#' Partial fishing mortality by fleet, gear or area
#'
#' Splits the annual fishing mortality among fleets or groups of fleets (gears, countries, areas), as the partial F
#' plots of the ICES / GFCM advice. Three methods:
#' \describe{
#'   \item{\code{"catch"}}{(default) the reported annual F (\code{ss3rep$exploitation$annual_F}, the
#'     \code{F_std} basis, e.g. Fbar over ages 1-3 as Z - M) split in proportion to the catch in numbers of each fleet
#'     over the F ages. Partial F's add up to the total F, also in area-structured models where fleets fish
#'     different parts of the stock.}
#'   \item{\code{"fbar"}}{the mean over the F ages of the annual F-at-age of each fleet (\code{\link{ss3_fa}()}, seasons
#'     weighted by duration), summed within groups. This is the F experienced by the fish in the area of the fleet:
#'     with \code{groups = "area"} it is the F by area, which does not add up to the stock F when the stock is spread
#'     across areas.}
#'   \item{\code{"apical"}}{the maximum over ages of the annual F-at-age of each fleet (as \code{apply(fa, 2, max)}),
#'     summed within groups.}
#' }
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param groups \code{NULL} (each catch fleet, default), \code{"area"} (fleets grouped by area) or a named list of
#'   fleet names or numbers, e.g. \code{list(Trawl = c("DK3a", "NOR3a"), Gillnet = "SWE3a")}
#' @param method \code{"catch"}, \code{"fbar"} or \code{"apical"} (see details)
#' @param ages ages over which F is averaged (default: the ages of the SS3 \code{F_std} basis; all ages if the basis
#'   is not age-based)
#' @param season keep the seasons (\code{"fbar"} and \code{"apical"} only; season F's on an annual scale)
#' @param years years to keep (default: model years)
#' @param long return the long data.table
#' @return FLQuants of partial F by group (units "f"), with the total F as attribute \code{"total"} (FLQuant) and the
#'   method as attribute \code{"method"}; or a data.table
#' @seealso \code{\link{ss3plot_partialF}()}, \code{\link{ss3_fa}()}, \code{\link{ss3_impact}()}
#' @examples
#' \dontrun{
#' pf <- ss3_partialF(ss3rep)                       # by fleet, adds up to F
#' ss3plot_partialF(pf)
#' ss3plot_partialF(ss3_partialF(ss3rep, groups = "area", method = "fbar"))   # F by area
#' }
#' @export
ss3_partialF <- function(ss3rep, groups = NULL, method = c("catch", "fbar", "apical"), ages = NULL,
                         season = FALSE, years = NULL, long = FALSE){
  method <- match.arg(method)
  fl <- ss3_fleets(ss3rep)
  fleets <- fl[type %in% c("catch", "bycatch")]$name
  grp <- .fleet_groups(ss3rep, groups, fleets)
  g <- data.table(name = unlist(grp), group = rep(names(grp), lengths(grp)))
  g[, fleet := match(name, .fleetnames(ss3rep))]
  yrs <- years %||% seq(ss3rep$startyr, ss3rep$endyr)
  ages <- ages %||% .fbar_ages(ss3rep)
  ex <- as.data.table(ss3rep$exploitation)
  tot <- if(nrow(ex) && "annual_F" %in% names(ex)) ex[Yr %in% yrs, .(data = annual_F[1]), by = .(year = Yr)] else NULL
  if(method == "catch"){
    if(season) warning("season = TRUE ignored for method = \"catch\"", call. = FALSE)
    if(is.null(tot)) stop("annual_F not found in ss3rep$exploitation", call. = FALSE)
    ca <- .catage_annual(ss3rep)[year %in% yrs]
    if(!is.null(ages)) ca <- ca[age %in% ages]
    ca <- merge(ca[, .(c = sum(data)), by = .(fleet, year)], g, by = "fleet")
    ca[, w := c/sum(c), by = year]
    d <- merge(ca[, .(w = sum(w)), by = .(group, year)], tot, by = "year")
    d <- d[, .(group, year, season = "all", data = w*data)]
  } else {
    fa <- ss3_fa(ss3rep, fleets = g$name, annual = !season, long = TRUE)[year %in% yrs]
    fa[, a := as.numeric(age)]
    fa <- fa[, .(data = mean(data)), by = .(name, year, season, a)]          # mean over sexes
    fa <- if(method == "fbar") fa[if(is.null(ages)) TRUE else a %in% ages, .(data = mean(data)), by = .(name, year, season)]
          else fa[, .(data = max(data)), by = .(name, year, season)]
    d <- merge(fa, g, by = "name")[, .(data = sum(data)), by = .(group, year, season)]
  }
  d <- d[, .(name = group, age = "all", year, unit = "unique", season, area = "unique", iter = "1", data)]
  if(long) return(d[])
  out <- long2flqs(d, quant = "age", units = "f", order = names(grp), years = yrs)
  if(!is.null(tot))
    attr(out, "total") <- long2flq(tot[, .(age = "all", year, unit = "unique", season = "all", area = "unique",
                                         iter = "1", data)], units = "f", years = yrs)
  attr(out, "method") <- method; attr(out, "ages") <- ages
  out
}

#' Plot partial F
#'
#' Partial fishing mortality by fleet, gear or area from \code{\link{ss3_partialF}()} as stacked areas (with the total
#' F as a line), lines, or shares.
#'
#' @param x FLQuants from \code{\link{ss3_partialF}()} (or any FLQuants of F by group)
#' @param type \code{"area"} (stacked, default for method "catch"), \code{"line"} (default otherwise) or
#'   \code{"share"}
#' @param total add the total F (attribute \code{"total"}) as a black line
#' @param ylab y-axis label
#' @return ggplot
#' @export
ss3plot_partialF <- function(x, type = NULL, total = TRUE, ylab = NULL){
  method <- attr(x, "method") %||% "catch"
  type <- match.arg(type %||% if(method == "catch") "area" else "line", c("area", "line", "share"))
  d <- .flq_long(x, 1, "F")
  bys <- length(unique(d$season)) > 1
  cols <- stats::setNames(ss3cols(nlevels(d$name)), levels(d$name))
  ylab <- ylab %||% switch(type, share = "Share of F",
                           paste0(switch(method, catch = "Partial F", fbar = "Fbar", apical = "Apical F"),
                                  if(!is.null(attr(x, "ages"))) paste0(" (ages ", min(attr(x, "ages")), "-",
                                                                         max(attr(x, "ages")), ")")))
  p <- ggplot2::ggplot(d, ggplot2::aes(year, data))
  p <- if(type == "line")
    p + ggplot2::geom_line(ggplot2::aes(colour = name), linewidth = 0.7) + ggplot2::scale_colour_manual(values = cols)
  else
    p + ggplot2::geom_area(ggplot2::aes(fill = name), position = if(type == "share") "fill" else "stack",
                           colour = "white", linewidth = 0.1, alpha = 0.9) + ggplot2::scale_fill_manual(values = cols)
  tf <- attr(x, "total")
  if(total && !is.null(tf) && type != "share" && !bys)
    p <- p + ggplot2::geom_line(data = .flq_long(FLCore::FLQuants(total = tf), 1, "F"), ggplot2::aes(year, data),
                                colour = "black", linewidth = 0.8, inherit.aes = FALSE) +
      ggplot2::labs(caption = "black line: total F (F_std)")
  if(bys) p <- p + ggplot2::facet_wrap(~ season, labeller = ggplot2::labeller(season = function(s) paste("season", s)))
  p + ggplot2::labs(x = "Year", y = ylab, fill = NULL, colour = NULL) + scale_x_year() + theme_ss3()
}

#' Plot F-at-age by fleet
#'
#' F-at-age curves of each fleet from \code{\link{ss3_fa}()}, averaged over the selected years, by season (colours)
#' and unit (line type). Shows the seasonal pattern and the selectivity of the fleets on the F scale.
#'
#' @param fa FLQuants from \code{\link{ss3_fa}()}
#' @param years years to average (default: the last three)
#' @param ncol number of panel columns
#' @return ggplot
#' @export
ss3plot_fa <- function(fa, years = NULL, ncol = NULL){
  d <- .flq_long(fa, 1, "F")
  years <- years %||% utils::tail(sort(unique(d$year)), 3)
  d <- d[year %in% years, .(data = mean(data)), by = .(name, quant, unit, season)]
  d[, age := as.numeric(as.character(quant))]
  # sex-invariant F-at-age: show one unit only
  if(length(unique(d$unit)) > 1){
    w <- dcast(d, name + age + season ~ unit, value.var = "data")
    if(all(vapply(w[, -(1:3)], function(v) isTRUE(all.equal(v, w[[4]])), TRUE))) d <- d[unit == unit[1]][, unit := "all"]
  }
  ns <- length(unique(d$season))
  cols <- if(ns > 1) stats::setNames(ss3seascols(ns), .lev(d$season)) else c(all = "black")
  d[, season := factor(season, levels = names(cols))]
  p <- ggplot2::ggplot(d, ggplot2::aes(age, data, colour = season, linetype = unit)) +
    ggplot2::geom_line(linewidth = 0.7) + ggplot2::scale_colour_manual(values = cols) +
    ggplot2::facet_wrap(~ name, ncol = ncol) +
    ggplot2::labs(x = "Age", y = "F at age (annual scale)", colour = "Season", linetype = "Unit",
                  subtitle = paste0("mean ", min(years), if(length(years) > 1) paste0("-", max(years)))) +
    theme_ss3() + ggplot2::theme(legend.title = ggplot2::element_text())
  if(length(unique(d$unit)) == 1) p <- p + ggplot2::guides(linetype = "none")
  p
}
