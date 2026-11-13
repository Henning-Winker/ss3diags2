# ensemble.R - model ensembles from several SS3 runs (as ss3diags::SSplotEnsemble)

# long table year, iter, SSB, Recr, Catch, F, stock (B/Btgt), harvest (F/Ftgt) of an FLStockR with iterations
.stk_draws <- function(stk){
  q <- ss3_flqs(stk, c("SSB", "Rec", "Catch", "F", "B/Btgt", "F/Ftgt"))
  v <- c(SSB = "SSB", Rec = "Recr", Catch = "Catch", F = "F", `B/Btgt` = "stock", `F/Ftgt` = "harvest")
  d <- Reduce(function(a, b) merge(a, b, by = c("year", "iter")), lapply(names(v), function(k){
    x <- as.data.table(methods::as(q[[k]], "data.frame"))[, .(year = as.numeric(as.character(year)),
                                                             iter = as.integer(as.character(iter)), data)]
    setnames(x, "data", v[[k]]); x }))
  d
}

#' Model ensemble from several SS3 runs
#'
#' Combines the MVLN draws of several SS3 runs into a model ensemble, as \code{ss3diags::SSdeltaMVLN(weight = )} with
#' \code{SSplotEnsemble(joint = TRUE)}: each model contributes draws in proportion to its weight, and the joint
#' ensemble pools them. Ratios to the targets (B/Btgt, F/Ftgt) are computed within each model, with its own
#' reference points, before pooling, so models with different benchmarks combine on the relative scale.
#'
#' @param x named list of SS3 reports (\code{r4ss::SS_output()}) or of FLStockR with iterations (e.g.
#'   \code{ss3stk(output = "iters")}), or FLStocks
#' @param weights model weights (default equal); normalised to sum to one
#' @param iters average number of draws per model: model \eqn{i} gets \code{round(iters * n * w_i)} draws for
#'   \eqn{n} models (FLStockR inputs are resampled to that number)
#' @param joint add the pooled ensemble (\code{"Joint"})
#' @param method,fwd,Fref passed to \code{\link{ss3stk}()} for SS3 reports
#' @param seed random seed
#' @param ... further arguments to \code{\link{ss3stk}()}
#' @return FLStocks: one FLStockR per model and the joint ensemble (iterations from all models; reference points with
#'   the iterations of the models; the ratios to the targets in \code{attr(, "status")}); the weights and draws per
#'   model as attribute \code{"weights"}
#' @seealso \code{\link{ss3plot_ensemble}()}
#' @examples
#' \dontrun{
#' ens <- ss3_ensemble(list(base = ss3rep1, highM = ss3rep2, lowh = ss3rep3), weights = c(0.5, 0.25, 0.25))
#' ss3plot_ensemble(ens)
#' }
#' @export
ss3_ensemble <- function(x, weights = NULL, iters = 1000, joint = TRUE, method = c("joint", "annual"), fwd = FALSE,
                         Fref = NULL, seed = 1, ...){
  method <- match.arg(method)
  if(is(x, "FLStock")) x <- list(x)
  x <- as.list(x)
  nm <- names(x); if(is.null(nm) || any(nm == "")) nm <- paste0("model", seq_along(x))
  names(x) <- nm
  w <- weights %||% rep(1, length(x))
  if(length(w) != length(x) || any(w < 0) || sum(w) <= 0) stop("'weights' must be non-negative, one per model", call. = FALSE)
  w <- w/sum(w); ni <- pmax(round(iters*length(x)*w), 0)
  old <- if(exists(".Random.seed", envir = globalenv())) get(".Random.seed", envir = globalenv()) else NULL
  on.exit(if(!is.null(old)) assign(".Random.seed", old, envir = globalenv()), add = TRUE)
  set.seed(seed)
  stks <- lapply(seq_along(x), function(i){
    if(ni[i] == 0) return(NULL)
    m <- x[[i]]
    if(.is_ss3rep(m))
      return(ss3stk(m, output = "iters", iters = ni[i], method = method, fwd = fwd, Fref = Fref, name = nm[i],
                    ...))
    if(!is(m, "FLStockR")) stop("models must be SS3 reports or FLStockR", call. = FALSE)
    nit <- FLCore::dims(m)$iter
    if(nit < 2) stop("FLStockR '", nm[i], "' has no iterations; use ss3stk(output = \"iters\")", call. = FALSE)
    id <- sample(seq_len(nit), ni[i], replace = ni[i] > nit)
    st <- ss3_flqs(m, c("B/Btgt", "F/Ftgt"))
    s <- FLCore::iter(m, id); FLCore::refpts(s) <- FLCore::iter(FLCore::refpts(m), if(FLCore::dims(FLCore::refpts(m))$iter > 1) id else 1)
    attributes(s)[c("basis", "endyr")] <- attributes(m)[c("basis", "endyr")]
    if(!is.null(attr(m, "status"))) attr(s, "status") <- list(stock = FLCore::iter(st[["B/Btgt"]], id),
                                                             harvest = FLCore::iter(st[["F/Ftgt"]], id))
    s
  })
  keep <- !vapply(stks, is.null, TRUE)
  stks <- stats::setNames(stks[keep], nm[keep])
  lab <- vapply(stks, function(s) paste(attr(s, "basis")$B %||% "Btgt", attr(s, "basis")$F %||% "Ftgt"), "")
  if(length(unique(lab)) > 1)
    warning("models have different targets (", paste(unique(lab), collapse = "; "),
            "): the joint ensemble combines B/Btgt and F/Ftgt on the relative scale", call. = FALSE)
  out <- stks
  if(joint && length(stks) > 1){
    yrs <- Reduce(intersect, lapply(stks, function(s) as.numeric(dimnames(s)$year)))
    off <- cumsum(c(0, vapply(stks, function(s) FLCore::dims(s)$iter, 1)))
    d <- rbindlist(lapply(seq_along(stks), function(i) .stk_draws(stks[[i]])[year %in% yrs][, iter := iter + off[i]]))
    js <- .stk_build(d, "Joint")
    rp <- lapply(seq_along(stks), function(i){
      r <- FLCore::refpts(stks[[i]]); n <- FLCore::dims(stks[[i]])$iter
      a <- array(NA_real_, dim = c(5, n), dimnames = list(params = c("Ftgt", "Btgt", "MSY", "B0", "R0"), iter = seq_len(n)))
      for(p in intersect(dimnames(a)$params, dimnames(r)$params)) a[p, ] <- rep_len(as.numeric(r[p, ]), n)
      a })
    ra <- do.call(cbind, rp); dimnames(ra) <- list(params = rownames(ra), iter = seq_len(ncol(ra)))
    js <- FLCore::FLStockR(js, refpts = FLCore::FLPar(ra))
    mk <- function(v) long2flq(d[, .(age = "1", year, iter, data = get(v))], quant = "age", units = "")
    attr(js, "status") <- list(stock = mk("stock"), harvest = mk("harvest"))
    attr(js, "basis") <- attr(stks[[1]], "basis")
    attr(js, "endyr") <- min(vapply(stks, function(s) attr(s, "endyr") %||% max(yrs), 1))
    out <- c(out, list(Joint = js))
  }
  out <- FLCore::FLStocks(out)
  attr(out, "weights") <- data.table(model = nm, weight = w, iters = ni)
  out
}

#' Plot a model ensemble
#'
#' Trajectories of a model ensemble from \code{\link{ss3_ensemble}()} with \code{\link{ss3plot_trj}()}: the models in
#' colour (median and quantile bands) and the joint ensemble in black, as \code{ss3diags::SSplotEnsemble()}.
#'
#' @param x FLStocks from \code{\link{ss3_ensemble}()}
#' @param type quantities, as in \code{\link{ss3plot_trj}()}
#' @param models show the individual models
#' @param joint show the joint ensemble
#' @param probs quantiles for the bands
#' @param colours model colours (default r4ss colours; the joint ensemble is black)
#' @param ... further arguments to \code{\link{ss3plot_trj}()}
#' @return ggplot
#' @export
ss3plot_ensemble <- function(x, type = c("SSB", "F", "Rec", "Catch", "B/Btgt", "F/Ftgt"), models = TRUE, joint = TRUE,
                             probs = c(0.05, 0.2, 0.5, 0.8, 0.95), colours = NULL, ...){
  nm <- names(x); mods <- setdiff(nm, "Joint")
  sel <- c(if(models) mods, if(joint && "Joint" %in% nm) "Joint")
  if(!length(sel)) stop("nothing to plot", call. = FALSE)
  cols <- c(stats::setNames(colours %||% ss3cols(length(mods)), mods), Joint = "black")[sel]
  p <- ss3plot_trj(FLCore::FLStocks(as.list(x)[sel]), type = type, probs = probs, colours = cols, ...)
  w <- attr(x, "weights")
  if(!is.null(w)) p <- p + ggplot2::labs(subtitle = paste0("weights: ", paste0(w$model, " ", round(w$weight, 2),
                                                                               collapse = ", ")))
  if(length(sel) == 1) p <- p + ggplot2::theme(legend.position = "none")
  p
}
