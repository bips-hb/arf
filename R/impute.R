
#' Missing value imputation with ARF
#' 
#' Imputed a dataset with missing values using adversarial random forests (ARF).
#' Calls \code{adversarial_rf}, \code{forde} and \code{expct}/\code{forge}.
#'
#' @param x Input data.
#' @param m Number of multiple imputations. The default is single imputation (\code{m=1}).
#' @param expectation Return expected value instead of multiple imputations. By default, for single imputation (\code{m=1}), the expected value is returned.
#' @param num_trees Number of trees in ARF.
#' @param min_node_size Minimum node size in ARF.
#' @param round Round imputed values to their respective maximum precision in the original data set?
#' @param finite_bounds Impose finite bounds on all continuous variables? See \code{\link{forde}}.
#' @param epsilon Slack parameter on empirical bounds; see \code{\link{forde}}.
#' @param verbose Print progress for \code{adversarial_rf}?
#' @param ... Extra parameters to be passed to \code{adversarial_rf}, \code{forde}
#'   and \code{expct}/\code{forge}.
#'
#' @return Imputed data. A single data table is returned for \code{m=1} and a list of data table for \code{m > 1}.
#' @export
#'
#' @examples
#' # Generate some missings
#' iris_na <- iris
#' for (j in 1:ncol(iris)) {
#'   iris_na[sample(1:nrow(iris), 5), j] <- NA
#' }
#' 
#' # Parallelization with doParallel
#' doParallel::registerDoParallel(cores = 2)
#' 
#' # Single imputation
#' iris_imputed <- arf::impute(iris_na, m = 1)
#' 
#' # Multiple imputation
#' iris_imputed <- arf::impute(iris_na, m = 20)
impute <- function(x, 
                   m = 1, 
                   expectation = ifelse(m==1, TRUE, FALSE), 
                   num_trees = 100L, 
                   min_node_size = 10L, 
                   round = TRUE, 
                   finite_bounds = "local",
                   epsilon = 1e-14, 
                   verbose = FALSE,
                   ...) {
  
  # To avoid data.table check issues
  idx <- . <- NULL
  
  if (m > 1 & expectation) {
    stop("Multiple imputation with expectation is not possible.")
  }
  if (!anyNA(x)) {
    message("No missing values found. Returning input data.")
    return(x)
  }
  if (utils::packageVersion("ranger") < "0.16.4") {
    stop("Imputation requires ranger >= 0.16.4. Consider updating the ranger package.")
  } 
  
  # Separate ... arguments for each function
  arg_names <- list(arf = names(as.list(args(adversarial_rf))), 
                    forde = names(as.list(args(forde))), 
                    forge = names(as.list(args(forge))),
                    expct = names(as.list(args(expct))))
  dot_args <- list(...)
  arf_args <- dot_args[names(dot_args) %in% arg_names$arf]
  forde_args <- dot_args[names(dot_args) %in% arg_names$forde]
  forge_args <- dot_args[names(dot_args) %in% arg_names$forge]
  expct_args <- dot_args[names(dot_args) %in% arg_names$expct]
  
  # ARF and FORDE
  arf <- do.call(adversarial_rf, c(x = list(x), 
                                   verbose = list(verbose), 
                                   num_trees = list(num_trees), 
                                   min_node_size = list(min_node_size),
                                   arf_args))
  psi <- do.call(forde, c(arf = list(arf), 
                          x = list(x), 
                          finite_bounds = list(finite_bounds), 
                          epsilon = list(epsilon),
                          forde_args))
  
  if (expectation) {
    # Expected value
    x_imputed <- do.call(expct, c(params = list(psi), 
                                  evidence = list(x), 
                                  round = list(round),
                                  expct_args))
  } else {
    # Multiple imputation
    x_synth <- do.call(forge, c(params = list(psi), 
                                n_synth = list(m), 
                                evidence = list(x), 
                                round = list(round),
                                forge_args))
    x_synth <- as.data.table(x_synth)
    x_synth[, idx := rep(1:m, nrow(x))]
    x_imputed <- split(x_synth, by = "idx")
    x_imputed <- lapply(x_imputed, function(x) x[, idx := NULL])
  }
  x_imputed
}