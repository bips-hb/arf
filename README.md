# arf: Adversarial random forests <a href='https://bips-hb.github.io/arf/'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/arf)](https://cran.r-project.org/package=arf)
[![check-standard](https://github.com/bips-hb/arf/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/bips-hb/arf/actions/workflows/check-standard.yaml)
<!-- badges: end -->

## Introduction
Adversarial random forests (ARFs) recursively partition data into fully factorized leaves, where features are jointly independent. The procedure is iterative, with alternating rounds of generation and discrimination. Data become increasingly realistic at each round, until original and synthetic samples can no longer be reliably distinguished. This is useful for several unsupervised learning tasks, such as density estimation and data synthesis. Methods for both are implemented in this package. ARFs naturally handle unstructured data with mixed continuous and categorical covariates. They inherit many of the benefits of RFs, including speed, flexibility, and solid performance with default parameters. 


## Installation
The `arf` package is available on `CRAN`:
```R
install.packages("arf")
```
To install the development version from GitHub using `devtools`, run:
```R
devtools::install_github("bips-hb/arf")
```

## Examples
Using Fisher's iris dataset, we train an ARF and estimate distribution parameters:
```R
# Train the ARF
arf <- adversarial_rf(iris)

# Estimate distribution parameters
psi <- forde(arf, iris)
```

### Density estimation
To estimate log-likelihoods:
```R
mean(lik(arf, psi, iris))
```

### Generative modeling
To generate 100 synthetic samples: 
```R
forge(psi, 100)
```

For more detailed examples, see the package [vignette](https://bips-hb.github.io/arf/articles/vignette.html).

## References
* Watson, D. S., Blesch, K., Kapar, J. & Wright, M. N. (2022). Adversarial random forests for density estimation and generative modeling. To appear in *Proceedings of the 26th International Conference on Artificial Intelligence and Statistics*. Preprint: https://arxiv.org/abs/2205.09435.
