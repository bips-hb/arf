pd <- "/home/burk/repos/bips-hb/arf"
scen <- commandArgs(TRUE)[1]
suppressMessages(pkgload::load_all(pd, quiet = TRUE))
seen <- character()
h <- function(expr) {
  withCallingHandlers(
    expr,
    message = function(m) {
      seen <<- c(seen, paste("MESSAGE:", trimws(conditionMessage(m))))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      seen <<- c(seen, paste("WARNING:", trimws(conditionMessage(w))))
      invokeRestart("muffleWarning")
    }
  )
}
pipeline <- function(par = TRUE) {
  arf <- h(adversarial_rf(iris, num_trees = 20, parallel = par, verbose = FALSE))
  psi <- h(forde(arf, iris, parallel = par))
  evi <- data.frame(Species = iris$Species[c(1, 51, 101, 2, 52, 102)])
  h(forge(psi, n_synth = 2, evidence = evi, stepsize = 2, parallel = par, verbose = FALSE))
  h(expct(psi, evidence = evi, stepsize = 2, parallel = par, verbose = FALSE))
  h(lik(psi, iris, arf = arf, batch = 50, parallel = par))
  h(forge(psi, n_synth = 5, parallel = par))
  invisible(NULL)
}
switch(
  scen,
  default = pipeline(),
  seq = pipeline(FALSE),
  verbose_off = {
    options(arf.verbose = FALSE)
    pipeline()
  },
  fork = {
    doParallel::registerDoParallel(cores = 2)
    pipeline()
  },
  mirai = {
    mirai::daemons(2)
    mirai::everywhere(suppressMessages(pkgload::load_all(pd, quiet = TRUE)), pd = pd)
    pipeline()
    mirai::daemons(0)
    pipeline()
  },
  mirai_nomori = {
    cat("mori installed:", requireNamespace("mori", quietly = TRUE), "\n")
    mirai::daemons(2)
    pipeline()
    mirai::daemons(0)
  },
  opt_mirai_no_daemons = {
    options(arf.backend = "mirai")
    r <- tryCatch(pipeline(), error = function(e) seen <<- c(seen, paste("ERROR:", conditionMessage(e))))
  },
  dofuture = {
    doFuture::registerDoFuture()
    future::plan(future::multisession, workers = 2)
    pipeline()
    future::plan("sequential")
  },
  future_mirai = {
    future::plan(future.mirai::mirai_multisession, workers = 2)
    mirai::everywhere(suppressMessages(pkgload::load_all(pd, quiet = TRUE)), pd = pd)
    pipeline()
    future::plan("sequential")
  },
  stop("unknown")
)
cat(sprintf("--- %s: %d condition(s)\n", scen, length(seen)))
for (s in unique(seen)) {
  cat("  ", sum(seen == s), "x ", s, "\n", sep = "")
}
