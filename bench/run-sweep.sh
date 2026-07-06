#!/usr/bin/env bash
# Kick off an arf backend sweep on a cluster node.
#
#   bench/run-sweep.sh                 # "clean": 1 thread per worker (fair scaling baseline)
#   bench/run-sweep.sh realistic       # "realistic": data.table + ranger at N threads each
#   bench/run-sweep.sh realistic 8     # ...with N = 8 threads
#
# By default runs bench/sweep.R (forde only). Set ARF_BENCH_SCRIPT to sweep the
# whole pipeline (forde/forge/expct/lik/adversarial_rf):
#   ARF_BENCH_SCRIPT=bench/sweep-ops.R ARF_BENCH_ITERS=3 bench/run-sweep.sh clean
#
# The two regimes answer different questions (see bench/README.md):
#   clean      -- N workers == N cores, so backend scaling is unconfounded.
#   realistic  -- data.table/ranger also multi-thread, as in naive real-world use;
#                 this intentionally oversubscribes at high worker counts, to test
#                 whether the clean interpretation still holds in practice.
#
# Any ARF_BENCH_* variable set in the environment overrides the defaults below.
set -euo pipefail
cd "$(dirname "$0")/.."   # repo root: sweep.R uses paths relative to it

# data.table's OpenMP threads busy-wait (spin) by default, which burns CPU while
# idle and inflates load average enormously under many workers -- and starves
# mirai's IPC coordination. Make idle threads sleep instead. Must be set before
# R starts; it propagates to the callr children and mirai daemons they spawn.
export OMP_WAIT_POLICY="${OMP_WAIT_POLICY:-passive}"

MODE="${1:-clean}"
NTHREADS="${2:-10}"

export ARF_BENCH_WORKERS="${ARF_BENCH_WORKERS:-1,2,4,8,16}"
export ARF_BENCH_N="${ARF_BENCH_N:-5000,20000}"
export ARF_BENCH_TREES="${ARF_BENCH_TREES:-100,200}"
export ARF_BENCH_ITERS="${ARF_BENCH_ITERS:-1}"

case "$MODE" in
  clean)
    export ARF_BENCH_DT_THREADS="${ARF_BENCH_DT_THREADS:-1}"
    export ARF_BENCH_RANGER_THREADS="${ARF_BENCH_RANGER_THREADS:-1}"
    ;;
  realistic)
    export ARF_BENCH_DT_THREADS="${ARF_BENCH_DT_THREADS:-$NTHREADS}"
    export ARF_BENCH_RANGER_THREADS="${ARF_BENCH_RANGER_THREADS:-$NTHREADS}"
    ;;
  *)
    echo "usage: $0 [clean|realistic] [nthreads]" >&2
    exit 1
    ;;
esac

SCRIPT="${ARF_BENCH_SCRIPT:-bench/sweep.R}"

echo "Mode: $MODE"
echo "  script         = $SCRIPT"
echo "  workers        = $ARF_BENCH_WORKERS"
echo "  n              = $ARF_BENCH_N"
echo "  trees          = $ARF_BENCH_TREES"
echo "  dt.threads     = $ARF_BENCH_DT_THREADS"
echo "  ranger.threads = $ARF_BENCH_RANGER_THREADS"
echo "  iters          = $ARF_BENCH_ITERS"
if [ "$SCRIPT" = "bench/sweep-ops.R" ]; then
  echo "  ops            = ${ARF_BENCH_OPS:-forde,forge,expct,lik,adversarial_rf (default)}"
  echo "  n_evidence     = ${ARF_BENCH_NEVIDENCE:-2000 (default)}"
  echo "  n_synth        = ${ARF_BENCH_NSYNTH:-1 (default)}"
  echo "  n_folds        = ${ARF_BENCH_NFOLDS:-8 (default)}"
fi
exec Rscript "$SCRIPT"
