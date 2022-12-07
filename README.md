<!-- badges: start -->
  [![check-standard](https://github.com/bips-hb/arf/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/bips-hb/arf/actions/workflows/check-standard.yaml)
  <!-- badges: end -->

## arf: Adversarial Random Forests
Adversarial random forests for density estimation and generative modeling

### Introduction
Adversarial random forests (ARFs) recursively partition data into fully factorized leaves, where features are jointly independent. The procedure is iterative, with alternating rounds of generation and discrimination. Data become increasingly realistic at each round, until original and synthetic samples can no longer be reliably distinguished. This is useful for several unsupervised learning tasks, such as density estimation and data synthesis. Methods for both are implemented in this package. ARFs naturally handle unstructured data with mixed continuous and categorical covariates. They inherit many of the benefits of RFs, including speed, flexibility, and solid performance with default parameters. 


### Installation
To install the development version from GitHub using `devtools`, run

```R
devtools::install_github("bips-hb/arf")
```

### Examples
#### Density estimation
To run an adversarial random forest on the iris data, perform density estimation and calculate log-likelihood on the same data:
```R
arf <- adversarial_rf(iris)
psi <- forde(arf, iris)
mean(lik(arf, psi, iris))
```

#### Generative modeling
To generate synthetic data based on the iris data: 
```R
arf <- adversarial_rf(iris)
psi <- forde(arf, iris)
forge(psi, 100)
```

### References
* Watson, D. S., Blesch, K., Kapar, J. & Wright, M. N. (2022). Adversarial random forests for density estimation and generative modeling. Preprint: https://arxiv.org/abs/2205.09435.
