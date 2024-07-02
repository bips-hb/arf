#' Shortcut likelihood function
#' 
#' Calls \code{adversarial_rf}, \code{forde} and \code{lik}.
#' For repeated application, it is faster to save outputs of \code{adversarial_rf}
#' and \code{forde} and pass them via \code{...} or directly use \code{lik}.
#' 
#' @param x Input data. Integer variables are recoded as ordered factors with
#'   a warning. See Details.
#' @param query Data frame of samples, optionally comprising just a subset of 
#'   training features. See Details of \code{lik}. Is set to \code{x} if \code{zero}. 
#' @param ... Extra parameters to be passed to \code{adversarial_rf}, \code{forde}
#'   and \code{lik}.
#'   
#' @return 
#' A vector of likelihoods, optionally on the log scale. A dataset of \code{n_synth} synthetic samples or of \code{nrow(x)} synthetic
#' samples if \code{n_synth} is undefined. 
#' 
#' 
#' @references 
#' Watson, D., Blesch, K., Kapar, J., & Wright, M. (2023). Adversarial random 
#' forests for density estimation and generative modeling. In \emph{Proceedings 
#' of the 26th International Conference on Artificial Intelligence and 
#' Statistics}, pp. 5357-5375.
#' 
#' 
#' @examples
#' # Estimate log-likelihoods
#' ll <- darf(iris)
#' 
#' # Partial evidence query
#' ll <- darf(iris, query = iris[1, 1:3])
#' 
#' # Condition on Species = "setosa"
#' ll <- darf(iris, query = iris[1, 1:3], evidence = data.frame(Species = "setosa"))
#' 
#' 
#' @seealso
#' \code{\link{arf}}, \code{\link{adversarial_rf}}, \code{\link{forde}}, \code{\link{forge}}
#' 
#'
#' @export
#' 

darf <- function(x, query = NULL, ...) {
  arg_names <- list(arf = names(as.list(args(adversarial_rf))), 
                    forde = names(as.list(args(forde))), 
                    lik = names(as.list(args(lik))))
  dot_args <- list(...)
  
  arf_args <- dot_args[names(dot_args) %in% arg_names$arf]
  forde_args <- dot_args[names(dot_args) %in% arg_names$forde]
  lik_args <- dot_args[names(dot_args) %in% arg_names$lik]
  
  if (!("verbose" %in% names(arf_args))) arf_args$verbose = F
  if (!("arf" %in% names(arf_args) | "params" %in% names(forde_args))) arf <- do.call(adversarial_rf, c(x = list(x), arf_args))
  
  if (!("params" %in% names(forde_args))) params <- do.call(forde, c(arf = list(arf), x = list(x), forde_args))
  
  if (is.null(query)) query <- x
  if (!("arf" %in% names(lik_args))) lik_args$arf <- arf
  do.call(lik, c(params = list(params),
                 query = list(query),
                 lik_args))

}


#' Shortcut sampling function
#' 
#' Calls \code{adversarial_rf}, \code{forde} and \code{forge}.
#' For repeated application, it is faster to save outputs of \code{adversarial_rf}
#' and \code{forde} and pass them via \code{...} or directly use \code{forge}.
#' 
#' @param x Input data. Integer variables are recoded as ordered factors with
#'   a warning. See Details.
#' @param n_synth Number of synthetic samples to generate. Is set to \code{nrow(x)} if
#' \code{NULL}.
#' @param ... Extra parameters to be passed to \code{adversarial_rf}, \code{forde}
#'   and \code{forge}.
#'   
#' @return 
#' A dataset of \code{n_synth} synthetic samples or of \code{nrow(x)} synthetic
#' samples if \code{n_synth} is undefined. 
#' 
#' 
#' @references 
#' Watson, D., Blesch, K., Kapar, J., & Wright, M. (2023). Adversarial random 
#' forests for density estimation and generative modeling. In \emph{Proceedings 
#' of the 26th International Conference on Artificial Intelligence and 
#' Statistics}, pp. 5357-5375.
#' 
#' 
#' @examples
#' # Generate 150 (size of original iris dataset) synthetic samples from the iris dataset
#' x_synth <- rarf(iris)
#' 
#' # Generate 100 synthetic samples from the iris dataset
#' x_synth <- rarf(iris, n_synth = 100)
#' 
#' # Condition on Species = "setosa"
#' x_synth <- rarf(iris, evidence = data.frame(Species = "setosa"))
#' 
#' @seealso
#' \code{\link{arf}}, \code{\link{adversarial_rf}}, \code{\link{forde}}, \code{\link{forge}}
#' 
#'
#' @export
#' 

rarf <- function(x, n_synth = NULL, ...) {
  arg_names <- list(arf = names(as.list(args(adversarial_rf))), 
                    forde = names(as.list(args(forde))), 
                    forge = names(as.list(args(forge))))
  dot_args <- list(...)
  
  arf_args <- dot_args[names(dot_args) %in% arg_names$arf]
  forde_args <- dot_args[names(dot_args) %in% arg_names$forde]
  forge_args <- dot_args[names(dot_args) %in% arg_names$forge]
  
  if (!("verbose" %in% names(arf_args))) arf_args$verbose = F
  if (!("arf" %in% names(arf_args) | "params" %in% names(forde_args))) arf <- do.call(adversarial_rf, c(x = list(x), arf_args))
  
  if (!("params" %in% names(forde_args))) params <- do.call(forde, c(arf = list(arf), x = list(x), forde_args))
  
  if(is.null(n_synth)) n_synth <- nrow(x)
  do.call(forge, c(params = list(params), 
                   n_synth = list(n_synth),
                   forge_args))
}


#' Shortcut expectation function
#' 
#' Calls \code{adversarial_rf}, \code{forde} and \code{expct}.
#' For repeated application, it is faster to save outputs of \code{adversarial_rf}
#' and \code{forde} and pass them via \code{...} or directly use \code{expct}.
#' 
#' @param x Input data. Integer variables are recoded as ordered factors with
#'   a warning. See Details.
#' @param ... Extra parameters to be passed to \code{adversarial_rf}, \code{forde}
#'   and \code{expct}.
#'   
#' @return 
#' A one row data frame with values for all query variables.
#' 
#' 
#' @references 
#' Watson, D., Blesch, K., Kapar, J., & Wright, M. (2023). Adversarial random 
#' forests for density estimation and generative modeling. In \emph{Proceedings 
#' of the 26th International Conference on Artificial Intelligence and 
#' Statistics}, pp. 5357-5375.
#' 
#' 
#' @examples
#' # What is the expected values of each feature?
#' earf(iris)
#' 
#' #' # What is the expected values of Sepal.Length?
#' earf(iris, query = "Sepal.Length")
#' 
#' # What if we condition on Species = "setosa"?
#' earf(iris, query = "Sepal.Length", evidence = data.frame(Species = "setosa"))
#' 
#' 
#' @seealso
#' \code{\link{arf}}, \code{\link{adversarial_rf}}, \code{\link{forde}}, \code{\link{expct}}
#' 
#'
#' @export
#' 

earf <- function(x, ...) {
  arg_names <- list(arf = names(as.list(args(adversarial_rf))), 
                    forde = names(as.list(args(forde))), 
                    expct = names(as.list(args(expct))))
  dot_args <- list(...)
  
  arf_args <- dot_args[names(dot_args) %in% arg_names$arf]
  forde_args <- dot_args[names(dot_args) %in% arg_names$forde]
  expct_args <- dot_args[names(dot_args) %in% arg_names$expct]
  
  if (!("verbose" %in% names(arf_args))) arf_args$verbose = F
  if (!("arf" %in% names(arf_args) | "params" %in% names(forde_args))) arf <- do.call(adversarial_rf, c(x = list(x), arf_args))
  
  if (!("params" %in% names(forde_args))) params <- do.call(forde, c(arf = list(arf), x = list(x), forde_args))
  
  do.call(expct, c(params = list(params),
                   expct_args))
  
}
