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
#'   \item{\code{arf.verbose}}{Report the selected backend once per backend
#'     configuration per session? Default \code{TRUE}; set \code{FALSE} to
#'     silence.}
#'   \item{\code{arf.block_rows}}{Cap on rows materialized per block of
#'     conditions in \code{\link{expct}}. Default \code{5e6}. Lower it to
#'     trade speed for memory on large forests with many conditions.}
#' }
#'
#' \code{arf.block_rows} does not affect results, only peak memory and speed.
#' For memory-constrained hardware, combine a small daemon count with a lower
#' \code{arf.block_rows} and the \code{batch}/\code{stepsize} arguments of
#' \code{\link{lik}}, \code{\link{forge}} and \code{\link{expct}}.
#'
#' The \code{"foreach"} backend uses whatever adapter is registered (e.g.
#' \code{doParallel}, \code{doFuture}). The \code{"mirai"} backend uses
#' \code{mirai} daemons and shares large read-only inputs (training data,
#' forest, learned parameters) across workers via \code{mori}, so workers do
#' not each copy them. Speed is comparable between the backends. The memory
#' benefit is largest for the tree-parallel operations (\code{\link{forde}},
#' \code{\link{adversarial_rf}}) on large forests with many workers, where
#' \code{foreach} memory grows with the worker count and \code{mirai} stays
#' much flatter (roughly half to a third at 16 workers in internal benchmarks).
#' For small workloads the daemon pool adds a fixed overhead that can outweigh
#' the sharing, so prefer \code{"mirai"} at scale and either backend otherwise.
#' Daemons started via \code{future.mirai} (e.g.
#' \code{plan(future.mirai::mirai_multisession)}) are detected like any other
#' mirai daemons, so futureverse users get the \code{"mirai"} backend
#' automatically. All backend packages are in Suggests; install the ones you
#' use.
#'
#' Reproducibility of stochastic operations (\code{\link{forge}}, categorical
#' \code{\link{expct}}) under parallel execution: \code{\link{set.seed}} only
#' governs the calling process, not the workers. With the \code{"mirai"}
#' backend, seed the daemons instead: \code{mirai::daemons(n, seed = 42)}
#' gives reproducible results provided the daemon count, the seed and the
#' sequence of calls on a fresh daemon pool are kept fixed (changing the
#' daemon count changes how work is chunked and therefore the random stream
#' assignment). This also applies to pools started via \code{future.mirai}:
#' \code{future}'s own seed machinery covers only work dispatched through its
#' API, which arf's backend bypasses, and \code{plan()} does not accept a
#' \code{seed} argument, so seed the pool with \code{mirai::daemons()}
#' directly. For the \code{"foreach"} backend, register \code{doRNG} on top of
#' the adapter (\code{doRNG::registerDoRNG(42)} after
#' \code{registerDoParallel()} or \code{registerDoFuture()}): results are then
#' reproducible, identical across adapters, and independent of the worker
#' count, provided \code{stepsize} is set explicitly (its default depends on
#' the worker count). Without \code{doRNG}, \code{doFuture} flags the
#' stochastic operations with "UNRELIABLE VALUE" warnings. Sequential
#' execution (\code{parallel = FALSE}) with \code{set.seed} is exact as
#' always.
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
#' # futureverse: future.mirai daemons are detected automatically
#' future::plan(future.mirai::mirai_multisession, workers = 4)
#' psi <- forde(arf, iris)
#'
#' # force a backend regardless of what is registered
#' options(arf.backend = "mirai")
#'
#' # silence the backend message
#' options(arf.verbose = FALSE)
#'
#' # reproducible parallel sampling with mirai: seeded daemons
#' # (fixed daemon count, fresh pool)
#' evi <- data.frame(Species = sample(levels(iris$Species), 100, replace = TRUE))
#' mirai::daemons(4, seed = 42)
#' # stepsize = evidence rows per step; 25 = 100 conditions / 4 workers,
#' # i.e. the default sizing, made explicit
#' x <- forge(psi, n_synth = 1, evidence = evi, stepsize = 25)
#' mirai::daemons(0)
#'
#' # reproducible parallel sampling with foreach: doRNG on top of the
#' # adapter; set stepsize explicitly (its default depends on worker count)
#' doParallel::registerDoParallel(cores = 4)
#' doRNG::registerDoRNG(42)
#' x <- forge(psi, n_synth = 1, evidence = evi, stepsize = 25)
#' }
#'
#' @name arf-options
#' @aliases arf.backend arf.verbose arf.block_rows
NULL
