#' arf package options
#'
#' Options controlling the parallel backend and its messaging, set via
#' \code{\link{options}}.
#'
#' @details
#' \describe{
#'   \item{\code{arf.backend}}{Parallel backend used when \code{parallel = TRUE}:
#'     \code{"foreach"} or \code{"mirai"}. If unset, arf uses \code{"mirai"} when
#'     mirai daemons are running and \code{"foreach"} otherwise.}
#'   \item{\code{arf.verbose}}{Report the selected backend once per session?
#'     Default \code{TRUE}; set \code{FALSE} to silence.}
#'   \item{\code{arf.chunk_factor}}{Chunks per worker for \code{mirai} tree
#'     maps (\code{\link{forde}}, \code{\link{adversarial_rf}} pruning).
#'     Default \code{1}. Values above 1 split the trees into smaller chunks,
#'     lowering each daemon's peak memory and smoothing uneven tree sizes, at
#'     the cost of more result transfers.}
#'   \item{\code{arf.block_rows}}{Cap on rows materialized per block of
#'     conditions in \code{\link{expct}}. Default \code{5e6}. Lower it to
#'     trade speed for memory on large forests with many conditions.}
#' }
#'
#' \code{arf.chunk_factor} and \code{arf.block_rows} do not affect results,
#' only peak memory and speed. For memory-constrained hardware, combine a
#' small daemon count with \code{arf.chunk_factor} above 1, a lower
#' \code{arf.block_rows}, and the \code{batch}/\code{stepsize} arguments of
#' \code{\link{lik}}, \code{\link{forge}} and \code{\link{expct}}.
#'
#' The \code{"foreach"} backend uses whatever adapter is registered (e.g.
#' \code{doParallel}, \code{doFuture}). The \code{"mirai"} backend uses
#' \code{mirai} daemons and shares the learned circuit across workers via
#' \code{mori}, so workers do not each copy it. Speed is comparable between the
#' backends. The memory benefit is largest for \code{\link{forde}} on large
#' forests with many workers, where \code{foreach} memory grows with the worker
#' count and \code{mirai} stays much flatter. For small workloads the daemon
#' pool adds a fixed overhead that can outweigh the sharing, so prefer
#' \code{"mirai"} at scale and either backend otherwise. All backend packages
#' are in Suggests; install the ones you use.
#'
#' @examples
#' \dontrun{
#' arf <- adversarial_rf(iris)
#'
#' # foreach backend
#' doParallel::registerDoParallel(cores = 4)
#' psi <- forde(arf, iris)
#'
#' # mirai backend: start daemons, then call as usual
#' mirai::daemons(4)
#' psi <- forde(arf, iris)
#' mirai::daemons(0)  # shut down when done
#'
#' # force a backend regardless of what is registered
#' options(arf.backend = "mirai")
#'
#' # silence the once-per-session backend message
#' options(arf.verbose = FALSE)
#' }
#'
#' @name arf-options
#' @aliases arf.backend arf.verbose
NULL
