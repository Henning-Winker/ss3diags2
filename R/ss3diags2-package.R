#' ss3diags2: Stock Synthesis diagnostics and data extraction in FLR
#'
#' Extraction functions read an \code{r4ss::SS_output()} list (\code{ss3rep}) and return FLR objects, or with
#' \code{long = TRUE} the underlying long \code{data.table} with the FLQuant dimension columns
#' (\code{len}/\code{age}, \code{year}, \code{unit}, \code{season}, \code{area}, \code{iter}, \code{data}) plus SS3
#' bookkeeping columns (\code{fleet}, \code{name}, ...).
#'
#' Observations:
#' \itemize{
#'   \item \code{\link{ss3_fleets}()}: fleet table (type, area, units, which data each fleet has)
#'   \item \code{\link{ss3_catch}()}: catch by fleet, \code{FLQuants}
#'   \item \code{\link{ss3_indices}()}: abundance indices, \code{FLIndices} of \code{FLIndexBiomass}
#'   \item \code{\link{ss3_lencomp}()}, \code{\link{ss3_agecomp}()}: compositions by fleet, \code{FLQuants}
#'   \item \code{\link{ss2jabba3}()}, \code{\link{ss2flicc}()}: the same, shaped as JABBA3 and FLicc inputs
#' }
#' @keywords internal
#' @import data.table
#' @importFrom methods is new as
#' @importFrom stats setNames weighted.mean
"_PACKAGE"
