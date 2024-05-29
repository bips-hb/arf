#' Quick-fix replacement function for vtruncnorm (because it seems to be buggy!)
#' 
#' @param a Vector of lower bounds of truncated Gaussians.
#' @param b Vector of higher bounds of truncated Gaussians.
#' @param mean Vector of means of corresponding untruncated Gaussians.
#' @param sd Vector of standard deviations of corresponding untruncated Gaussians.
#' 
#' @return Vector of variances of truncated Gaussians.
#' 
#' @import data.table
#' @importFrom truncnorm vtruncnorm rtruncnorm
#' @importFrom stats dnorm pnorm
#' @importFrom stats var sd
#' @keywords internal

var_truncnorm <- function(a = -Inf, b = Inf, mean = 0, sd = 1) {
  
  var <- NULL
  
  var_t <- vtruncnorm(a,b,mean,sd)
  if (any(is.na(var_t) | is.infinite(var_t))) {
    idx_rpl <- which(is.na(var_t) | is.infinite(var_t))
    a <- a[idx_rpl]
    b <- b[idx_rpl]
    mean <- mean[idx_rpl]
    sd <- sd[idx_rpl]
    alpha <- (a - mean)/sd
    beta <- (b - mean)/sd
    Z <- pnorm(beta) - pnorm(alpha)
    dn_alpha <- dnorm(alpha)
    dn_beta <- dnorm(beta)
    dn_prod_alpha <- fifelse(is.finite(alpha), alpha*dn_alpha, 0)
    dn_prod_beta <- fifelse(is.finite(beta), beta*dn_beta, 0)
    var_t[idx_rpl] <- sd^2*(1- (dn_prod_beta - dn_prod_alpha)/Z - ((dn_alpha - dn_beta)/Z)^2)
    if(any(is.na(var_t))) {
      idx_rpl <- which(is.na(var_t))
      a <- a[idx_rpl]
      b <- b[idx_rpl]
      mean <- mean[idx_rpl]
      sd <- sd[idx_rpl]
      n <- 1e6
      var_sample_dt <- data.table(a = a, b = b, mean = mean, sd = sd)
      var_sample_dt[, var := apply(.SD, 1, \(x) var(rtruncnorm(1e6, x["a"], x["b"], x["mean"] ,x["sd"]))*(n-1)/n)]
      var_t[idx_rpl] <- var_sample_dt[, var]
    }
  }
  var_t
}

#' Function for calculating truncated normal parameters for given bounds and 
#' desired mean and standard deviation.
#' 
#' @param a Vector of lower bounds of truncated Gaussians.
#' @param b Vector of higher bounds of truncated Gaussians.
#' @param mean Vector of desired means of truncated Gaussians.
#' @param sd Vector of desired standard deviations of truncated Gaussians.
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel}.
#' 
#' @return Parameter values for location and scale parameters for desired conditions.
#' 
#' @import data.table
#' @importFrom foreach foreach %do% %dopar% 
#' @importFrom truncnorm etruncnorm
#' @importFrom stats optim
#' @keywords internal

truncnorm_params <- function(a, b, mean, sd, parallel = TRUE) {
  
  i <- NULL
  
  doOptim <- function(a, b, mean, sd){ 
    init <- c(mean, log(sd*2))
    
    f <- function(params) {
      mu <- params[1]
      sigma <- exp(params[2])
      vt <- var_truncnorm(a, b, mu, sigma)
      if (vt > 0) {
        dev_sd <- (sd - sqrt(vt))
      } else {
        dev_sd <- sd^2 - vt
      }
      (mean - etruncnorm(a,b,mu,sigma))^2 + dev_sd^2
    }
    # solve for mu and sigma
    params <- optim(init, f)$par
    data.table(mu = params[1], sigma = (exp(params[2])))
  }
  len_a <- length(a)
  if (len_a > 1) {
    if (parallel) {
      result <- foreach(i = 1:len_a, .combine = rbind) %dopar% doOptim(a[i],b[i],mean[i],sd[i])
    } else {
      result <- foreach(i = 1:len_a, .combine = rbind) %do% doOptim(a[i],b[i],mean[i],sd[i])
    }
  } else result <- doOptim(a,b,mean,sd)
  
  result
  
}

#' Function for performing MLE for truncated Gaussians.
#' 
#' @param x Vector of input data.
#' @param a Vector of lower bounds of truncated Gaussians.
#' @param b Vector of higher bounds of truncated Gaussians.
#' 
#' @return Parameter values for location and scale parameters for desired conditions.
#' 
#' @import data.table
#' @importFrom truncnorm etruncnorm
#' @importFrom stats optim
#' @keywords internal

truncnorm_params_mle <- function(x, a = -Inf, b = Inf) {
  
  init <- c(mean(x), sd(x))
  
  f <- function(params) {
    mu <- params[1]
    sigma <- (params[2])
    
    -sum(log(dtruncnorm(x, a,b,mu,sigma)))
  }
  # solve for mu and sigma
  params <- optim(init, f)$par
  data.table(mu = params[1], sigma = params[2])
}

#' Function for correcting forde output with the function truncnorm_params.
#' 
#' @param params Circuit parameters learned via \code{\link{forde}}.
#' @param threshold Threshold value for deviations of mean and standard deviation
#'   calculated from input parameters under which no correction will be performed.
#' @param discard_deteriorations Discard correction where sum of the deviation of
#' mean and standard deviation is higher than with input parameters?
#' @param keep_old_params Include input parameters in output?
#' @param show_deviaton Include calculated deviations in output?
#' @param parallel Compute in parallel? Must register backend beforehand, e.g. 
#'   via \code{doParallel}.
#' 
#' @return Input circuit parameters with corrected mu and sigma.
#' 
#' @examples
#' # Train ARF and estimate leaf parameters
#' arf <- adversarial_rf(iris)
#' params <- forde(arf, iris)
#' 
#' # Correct params$cnt values for mu and sigma
#' params_corrected <- arf:::correct_params(params)
#' 
#' @import data.table
#' @importFrom truncnorm etruncnorm
#' @keywords internal

correct_params <- function(params, threshold = 1e-5, discard_deteriorations = T, keep_old_params = F, show_deviation = F, parallel = T) {
  
  . <- mu <- sigma <- mu_old <- sigma_old <- dev_mu_old <- dev_sigma_old <- dev_mu <- dev_sigma <- NULL
  
  params_correct <- copy(params)
  cnt <- copy(params_correct$cnt)
  cnt[, `:=` (mu_old = mu, sigma_old = sigma)]
  cnt[, c("dev_mu_old", "dev_sigma_old") := .(etruncnorm(min,max,mu_old, sigma_old) - mu_old, sqrt(var_truncnorm(min,max, mu_old, sigma_old)) -sigma_old)]
  cnt[abs(dev_mu_old) > threshold | abs(dev_sigma_old) > threshold, c("mu","sigma") := truncnorm_params(min, max, mean = mu, sd = sigma, parallel = parallel)]
  if (show_deviation | discard_deteriorations) {
    cnt[, c("dev_mu","dev_sigma") := .(etruncnorm(min,max,mu, sigma) - mu_old, sqrt(var_truncnorm(min,max, mu, sigma)) -sigma_old)]
    if (discard_deteriorations) {
      cnt[abs(dev_mu) + abs(dev_sigma) >= abs(dev_mu_old) + abs(dev_sigma_old), c("mu", "sigma") := .(mu_old, sigma_old)]
    }
  } else {
    cnt[, c("dev_mu_old", "dev_sigma_old") := NULL]
  }
  if (!keep_old_params) {
    cnt[, `:=` (mu_old = NULL, sigma_old = NULL)]
  }
  params_correct$cnt <- cnt[]
  params_correct
}