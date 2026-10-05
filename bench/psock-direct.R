suppressMessages(pkgload::load_all(
  "/tmp/claude-1000/-home-burk-repos-bips-hb-arf/32a9277e-c06c-4514-b8a7-44e19959cf92/scratchpad/mori-fe-wt",
  quiet = TRUE
))
options(arf.verbose = FALSE)
data.table::setDTthreads(1)
set.seed(1)
n <- 10000
p <- 30
X <- data.frame(matrix(rnorm(n * p), n))
X$grp <- factor(sample(letters[1:5], n, TRUE))
set.seed(2)
arf <- adversarial_rf(X, num_trees = 100, parallel = FALSE, verbose = FALSE)
hwm <- function(pid) {
  l <- readLines(sprintf("/proc/%d/status", pid))
  as.numeric(sub("VmHWM:\\s+(\\d+) kB", "\\1", grep("VmHWM", l, value = TRUE))) / 1024
}
shares <- 0
for (variant in c("plain", "mori")) {
  Sys.setenv(ARF_EXP_MORI_FOREACH = if (variant == "mori") "1" else "")
  cl <- parallel::makeCluster(4)
  parallel::clusterCall(
    cl,
    function(p) {
      suppressMessages(pkgload::load_all(p, quiet = TRUE))
      data.table::setDTthreads(1)
    },
    "/tmp/claude-1000/-home-burk-repos-bips-hb-arf/32a9277e-c06c-4514-b8a7-44e19959cf92/scratchpad/mori-fe-wt"
  )
  doParallel::registerDoParallel(cl)
  pids <- unlist(parallel::clusterCall(cl, Sys.getpid))
  n_shm0 <- length(list.files("/dev/shm", all.files = TRUE, no.. = TRUE))
  t <- system.time(psi <- forde(arf, X, parallel = TRUE))[["elapsed"]]
  cat(sprintf(
    "%-5s share() calls=%d  time=%.1fs  worker peak RSS (MB): %s  sum=%.0f  parent peak=%.0f\n",
    variant,
    shares,
    t,
    paste(round(sapply(pids, hwm)), collapse = "/"),
    sum(sapply(pids, hwm)),
    hwm(Sys.getpid())
  ))
  shares <- 0
  parallel::stopCluster(cl)
}
