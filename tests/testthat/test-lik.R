# Independent brute-force mixture density: sum over all leaves of
# wt * prod_j p_j(x_j), with no short-cuts and no filtering, so a leaf where one
# variable has zero density contributes zero. Deliberately not sharing code with
# lik(). truncnorm only.
lik_ref <- function(psi, query, log = TRUE) {
  query <- as.data.frame(query)
  f_idx <- as.character(psi$forest$f_idx)
  wt <- psi$forest$cvg / max(psi$forest$tree)
  out <- vapply(
    seq_len(nrow(query)),
    function(i) {
      dens <- stats::setNames(rep(1, length(f_idx)), f_idx)
      for (v in colnames(query)) {
        xv <- query[i, v]
        p <- stats::setNames(rep(0, length(f_idx)), f_idx)
        if (is.factor(xv)) {
          prob <- psi$cat[variable == v & val == as.character(xv)]
          p[as.character(prob$f_idx)] <- prob$prob
        } else {
          cnt <- psi$cnt[variable == v]
          d <- truncnorm::dtruncnorm(xv, a = cnt$min, b = cnt$max, mean = cnt$mu, sd = cnt$sigma)
          d[xv == cnt$min] <- 0
          p[as.character(cnt$f_idx)] <- d
        }
        dens <- dens * p
      }
      sum(wt * dens)
    },
    numeric(1)
  )
  if (isTRUE(log)) log(out) else out
}

make_mixed <- function(n) {
  data.frame(
    matrix(stats::rnorm(n * 2), n),
    g = factor(sample(letters[1:5], n, TRUE), levels = letters[1:5]),
    h = factor(sample(c("u", "v"), n, TRUE), levels = c("u", "v"))
  )
}

set.seed(1)
dat <- make_mixed(300)
new_dat <- make_mixed(20)
arf <- adversarial_rf(dat, num_trees = 5, parallel = FALSE, verbose = FALSE)
psi <- forde(arf, dat, parallel = FALSE)

test_that("mixed partial evidence matches brute-force mixture", {
  q <- dat[1:5, c("X1", "g")]
  expect_equal(lik(psi, q, parallel = FALSE), lik_ref(psi, q))
})

test_that("mixed total evidence matches brute-force mixture, with and without arf", {
  x <- new_dat[1:10, ]
  expect_equal(lik(psi, x, arf = arf, parallel = FALSE), lik_ref(psi, x))
  expect_equal(suppressWarnings(lik(psi, x, parallel = FALSE)), lik_ref(psi, x))
})

test_that("mixed likelihoods do not depend on the rest of the batch", {
  q <- dat[1:5, c("X1", "X2", "g")]
  per_row <- vapply(seq_len(nrow(q)), function(i) lik(psi, q[i, ], parallel = FALSE), numeric(1))
  expect_equal(lik(psi, q, parallel = FALSE), per_row)
  expect_equal(lik(psi, q, batch = 2, parallel = FALSE), per_row)
  expect_equal(lik(psi, q[5:1, ], parallel = FALSE), rev(per_row))
})

test_that("duplicate categorical patterns get their own likelihoods", {
  # Non-adjacent duplicates: row 3 repeats row 1's factor pattern
  q <- data.frame(
    X1 = dat$X1[1:4],
    X2 = dat$X2[1:4],
    g = factor(c("a", "b", "a", "c"), levels = letters[1:5]),
    h = factor(c("u", "v", "u", "v"), levels = c("u", "v"))
  )
  expect_equal(suppressWarnings(lik(psi, q, parallel = FALSE)), lik_ref(psi, q))

  q_cat <- q[, c("g", "h")]
  expect_equal(lik(psi, q_cat, parallel = FALSE), lik_ref(psi, q_cat))
})

test_that("pure continuous and pure categorical queries match brute-force mixture", {
  q_cnt <- new_dat[1:10, c("X1", "X2")]
  expect_equal(lik(psi, q_cnt, parallel = FALSE), lik_ref(psi, q_cnt))

  q_cat <- new_dat[1:10, c("g", "h")]
  expect_equal(lik(psi, q_cat, parallel = FALSE), lik_ref(psi, q_cat))
})

test_that("impossible mixed rows get zero likelihood", {
  q <- dat[1:3, c("X1", "g")]
  q$X1[2] <- 1e6
  expect_equal(lik(psi, q, parallel = FALSE), lik_ref(psi, q))
  expect_equal(lik(psi, q, parallel = FALSE)[2], -Inf)
})
