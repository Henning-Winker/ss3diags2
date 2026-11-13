# francis.R - Francis (2011) TA1.8 data weights for composition data

.ta18_comps <- function(ss3rep, quant, fleets, part, seas){
  db <- .comp_db(ss3rep, quant)
  nm <- .fleetnames(ss3rep)
  fl <- .pick_fleets(ss3rep, fleets, nm[sort(unique(db$Fleet))])
  db <- db[nm[Fleet] %in% fl & .col(db, "Part", 0) %in% part]
  if(!nrow(db)) return(NULL)
  db[, `:=`(sexes = .col(db, "Sexes", 0), lbin = .col(db, "Lbin_lo", NA_real_), method = .col(db, "method", 1),
            sea = if(seas == "sep") Seas else 0L)]
  # one composition: fleet x year x season (if separate) x sex option x method (size) x length bin (con)
  by <- c("Fleet", "Yr", "sea", if(quant == "con") "lbin", if(quant != "con") "sexes", if(quant == "size") "method")
  db[, .(obs = sum(Obs*Bin)/sum(Obs), exp = sum(Exp*Bin)/sum(Exp),
         se = sqrt((sum(Exp*Bin^2)/sum(Exp) - (sum(Exp*Bin)/sum(Exp))^2)/mean(Nsamp_adj)),
         n = mean(Nsamp_adj), seas = mean(Seas)), by = by]
}

#' Francis TA1.8 weights for composition data
#'
#' Data weights for length, age, conditional age-at-length or generalised size compositions after method TA1.8 of
#' Francis (2011), as \code{r4ss::SSMethod.TA1.8()} and \code{ss3diags::SScompsTA1.8()}: the mean length (or age) of
#' each composition is compared with its expected value, the standardised residuals use the variance implied by the
#' adjusted sample size, and the weight is \eqn{w = 1/\mathrm{var}(\mathrm{std.res})}. \eqn{w < 1} means the
#' compositions are overweighted (the mean lengths vary more than the sample sizes imply): multiply the SS3 variance
#' adjustment of the fleet by \eqn{w}. The 95\% interval is a bootstrap of the residuals.
#'
#' @param ss3rep output of \code{r4ss::SS_output()}
#' @param quant \code{"len"} (default), \code{"age"}, \code{"con"} (conditional age-at-length: one composition per
#'   length bin) or \code{"size"}
#' @param fleets fleets (numbers or names; default all with data)
#' @param part partitions combined (0 whole, 1 discard, 2 retained; default all)
#' @param seas \code{"comb"} (seasons combined in a composition, default) or \code{"sep"}
#' @param boot bootstrap samples for the interval
#' @param seed random seed for the bootstrap
#' @return data.table by fleet: \code{w} (Francis weight), \code{lo}, \code{hi} (95\% bootstrap interval), \code{n}
#'   (number of compositions), \code{var_adj} (current SS3 variance adjustment, from the comp fit summary of the
#'   report) and \code{var_adj.new} (\code{var_adj * w}, the input for the next run; not meaningful for
#'   Dirichlet-multinomial fleets, whose weight is estimated); the mean lengths / ages with their intervals as attribute \code{"comps"} (for
#'   \code{\link{ss3plot_francis}()})
#' @references Francis, R.I.C.C. (2011). Data weighting in statistical fisheries stock assessment models. Canadian
#'   Journal of Fisheries and Aquatic Sciences 68: 1124-1138.
#' @examples
#' data(ss3simple)
#' ss3_francis(ss3simple)
#' ss3_francis(ss3simple, "age")
#' @export
ss3_francis <- function(ss3rep, quant = c("len", "age", "con", "size"), fleets = NULL, part = 0:2,
                        seas = c("comb", "sep"), boot = 1000, seed = 1){
  quant <- match.arg(quant); seas <- match.arg(seas)
  pl <- .ta18_comps(ss3rep, quant, fleets, part, seas)
  if(is.null(pl)) stop("no compositions for these fleets / partitions", call. = FALSE)
  pl[, std.res := (obs - exp)/se]
  old <- if(exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if(!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed)
  nm <- .fleetnames(ss3rep)
  fs <- c(len = "Length_Comp_Fit_Summary", age = "Age_Comp_Fit_Summary", con = "Age_Comp_Fit_Summary",
          size = "Size_Comp_Fit_Summary")[[quant]]
  va <- as.data.table(ss3rep[[fs]] %||% data.frame())
  out <- rbindlist(lapply(sort(unique(pl$Fleet)), function(f){
    r <- pl[Fleet == f & is.finite(std.res)]$std.res
    if(length(r) < 2) return(NULL)
    w <- 1/stats::var(r)
    bt <- vapply(seq_len(boot), function(i) 1/stats::var(sample(r, replace = TRUE)), 1)
    ci <- stats::quantile(bt[is.finite(bt)], c(0.025, 0.975), names = FALSE)
    vadj <- NA_real_
    if(nrow(va) && all(c("Fleet", "Curr_Var_Adj") %in% names(va))){
      v <- suppressWarnings(as.numeric(va[as.numeric(Fleet) == f]$Curr_Var_Adj)); if(length(v)) vadj <- v[1]
    }
    data.table(fleet = nm[f], w = w, lo = ci[1], hi = ci[2], n = length(r), var_adj = vadj, var_adj.new = vadj*w)
  }))
  pl[, `:=`(name = nm[Fleet], w = out$w[match(nm[Fleet], out$fleet)])]
  pl[, `:=`(lo = obs - 2*se, hi = obs + 2*se, lo.adj = obs - 2*se/sqrt(w), hi.adj = obs + 2*se/sqrt(w))]
  setattr(out, "comps", pl[]); setattr(out, "quant", quant)
  out[]
}

#' Plot Francis TA1.8 mean lengths or ages
#'
#' Observed mean length (or age) of each composition with \eqn{\pm 2} SE from the adjusted sample size (thick bars)
#' and after the Francis weighting (thin bars), and the expected mean (line), by fleet over years (by year over length
#' bins for conditional age-at-length), with the Francis weight and its interval in the panel title.
#'
#' @param x output of \code{\link{ss3_francis}()}
#' @param ncol number of panel columns
#' @return ggplot
#' @export
ss3plot_francis <- function(x, ncol = NULL){
  pl <- attr(x, "comps"); quant <- attr(x, "quant") %||% "len"
  if(is.null(pl)) stop("x must be the output of ss3_francis()", call. = FALSE)
  pl <- copy(pl)
  lab <- x[, .(name = fleet, lab = sprintf("%s: w = %.3f (%.3f-%.3f)", fleet, w, lo, hi))]
  pl <- merge(pl, lab, by = "name")
  con <- quant == "con"
  pl[, `:=`(xx = if(con) lbin else Yr, panel = if(con) paste(lab, Yr) else lab)]
  grp <- intersect(c("sexes", "method"), names(pl))
  pl[, grp := do.call(paste, c(list(panel), .SD)), .SDcols = grp]
  ggplot2::ggplot(pl, ggplot2::aes(xx)) +
    ggplot2::geom_linerange(ggplot2::aes(ymin = lo, ymax = hi), linewidth = 1.2, colour = "grey45") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo.adj, ymax = hi.adj), width = 0.3, linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(y = obs), shape = 21, fill = "grey85", size = 1.8) +
    ggplot2::geom_line(ggplot2::aes(y = exp, group = grp), colour = "#0000CC", linewidth = 0.6) +
    ggplot2::facet_wrap(~ panel, scales = "free", ncol = ncol) +
    ggplot2::labs(x = if(con) "Length bin" else "Year",
                  y = paste("Mean", switch(quant, len = "length", size = "size", "age")),
                  caption = "thick: +/- 2 SE (adjusted N); thin: after Francis weighting; line: expected") +
    theme_ss3()
}
