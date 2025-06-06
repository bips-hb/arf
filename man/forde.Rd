% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forde.R
\name{forde}
\alias{forde}
\title{Forests for Density Estimation}
\usage{
forde(
  arf,
  x,
  oob = FALSE,
  family = "truncnorm",
  finite_bounds = c("no", "local", "global"),
  alpha = 0,
  epsilon = 0,
  parallel = TRUE
)
}
\arguments{
\item{arf}{Pre-trained \code{\link{adversarial_rf}}. Alternatively, any
object of class \code{ranger}.}

\item{x}{Training data for estimating parameters.}

\item{oob}{Only use out-of-bag samples for parameter estimation? If
\code{TRUE}, \code{x} must be the same dataset used to train \code{arf}.
Set to \code{"inbag"} to only use in-bag samples. Default is \code{FALSE},
i.e. use all observations.}

\item{family}{Distribution to use for density estimation of continuous
features. Current options include truncated normal (the default
\code{family = "truncnorm"}) and uniform (\code{family = "unif"}). See
Details.}

\item{finite_bounds}{Impose finite bounds on all continuous variables? If
\code{"local"}, infinite bounds are set to empirical extrema within leaves.
If \code{"global"}, infinite bounds are set to global empirical extrema.
if \code{"no"} (the default), infinite bounds are left unchanged.}

\item{alpha}{Optional pseudocount for Laplace smoothing of categorical
features. This avoids zero-mass points when test data fall outside the
support of training data. Effectively parameterizes a flat Dirichlet prior
on multinomial likelihoods.}

\item{epsilon}{Optional slack parameter on empirical bounds when
\code{finite_bounds != "no"}. This avoids zero-density points when test
data fall outside the support of training data. The gap between lower and
upper bounds is expanded by a factor of \code{1 + epsilon}.}

\item{parallel}{Compute in parallel? Must register backend beforehand, e.g.
via \code{doParallel} or \code{doFuture}; see examples.}
}
\value{
A \code{list} with 5 elements: (1) parameters for continuous data; (2)
parameters for discrete data; (3) leaf indices and coverage; (4) metadata on
variables; and (5) the data input class. This list is used for estimating
likelihoods with \code{\link{lik}} and generating data with \code{\link{forge}}.
}
\description{
Uses a pre-trained ARF model to estimate leaf and distribution parameters.
}
\details{
\code{forde} extracts leaf parameters from a pretrained forest and learns
distribution parameters for data within each leaf. The former includes
coverage (proportion of data falling into the leaf) and split criteria. The
latter includes proportions for categorical features and mean/variance for
continuous features. The result is a probabilistic circuit, stored as a
\code{data.table}, which can be used for various downstream inference tasks.

Currently, \code{forde} only provides support for a limited number of
distributional families: truncated normal or uniform for continuous data,
and multinomial for discrete data.

Though \code{forde} was designed to take an adversarial random forest as
input, the function's first argument can in principle be any object of class
\code{ranger}. This allows users to test performance with alternative
pipelines (e.g., with supervised forest input). There is also no requirement
that \code{x} be the data used to fit \code{arf}, unless \code{oob = TRUE}.
In fact, using another dataset here may protect against overfitting. This
connects with Wager & Athey's (2018) notion of "honest trees".
}
\examples{
# Train ARF and estimate leaf parameters
arf <- adversarial_rf(iris)
psi <- forde(arf, iris)

# Generate 100 synthetic samples from the iris dataset
x_synth <- forge(psi, n_synth = 100)

# Condition on Species = "setosa" and Sepal.Length > 6
evi <- data.frame(Species = "setosa",
                  Sepal.Length = "(6, Inf)")
x_synth <- forge(psi, n_synth = 100, evidence = evi)

# Estimate average log-likelihood
ll <- lik(psi, iris, arf = arf, log = TRUE)
mean(ll)

# Expectation of Sepal.Length for class setosa
evi <- data.frame(Species = "setosa")
expct(psi, query = "Sepal.Length", evidence = evi)

\dontrun{
# Parallelization with doParallel
doParallel::registerDoParallel(cores = 4)

# ... or with doFuture
doFuture::registerDoFuture()
future::plan("multisession", workers = 4)
}


}
\references{
Watson, D., Blesch, K., Kapar, J., & Wright, M. (2023). Adversarial random
forests for density estimation and generative modeling. In \emph{Proceedings
of the 26th International Conference on Artificial Intelligence and
Statistics}, pp. 5357-5375.

Wager, S. & Athey, S. (2018). Estimation and inference of heterogeneous
treatment effects using random forests. \emph{J. Am. Stat. Assoc.},
\emph{113}(523): 1228-1242.
}
\seealso{
\code{\link{arf}}, \code{\link{adversarial_rf}}, \code{\link{forge}},
\code{\link{expct}}, \code{\link{lik}}
}
