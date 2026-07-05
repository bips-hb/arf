#!/usr/bin/env bash
# Kick off the forde() backend sweep on a cluster node.
#
#   bench/run-sweep.sh                 # "clean": 1 thread per worker (fair scaling baseline)
#   bench/run-sweep.sh realistic       # "realistic": data.table + ranger at N threads each
#   bench/run-sweep.sh realistic 8     # ...with N = 8 threads
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

echo "Mode: $MODE"
echo "  workers        = $ARF_BENCH_WORKERS"
echo "  n              = $ARF_BENCH_N"
echo "  trees          = $ARF_BENCH_TREES"
echo "  dt.threads     = $ARF_BENCH_DT_THREADS"
echo "  ranger.threads = $ARF_BENCH_RANGER_THREADS"
echo "  iters          = $ARF_BENCH_ITERS"
exec Rscript bench/sweep.R
