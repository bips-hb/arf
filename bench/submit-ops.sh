#!/usr/bin/env bash
# Submit one slurm job per benchmark config so ops run in parallel on the
# cluster. Each job runs bench/run-sweep.sh clean with op-specific parameters;
# forge/expct get small/large and separate/or variants ("or" delegates the
# parallelism to cforde, so those configs benchmark the cforde backend).
#
#   bench/submit-ops.sh            # submit all configs
#   DRY_RUN=1 bench/submit-ops.sh  # print sbatch commands without submitting
#
# Slurm sizing (env-overridable):
#   ARF_BENCH_SLURM_CPUS  cpus per job  (default 17: max workers 16 + orchestrator)
#   ARF_BENCH_SLURM_MEM   memory        (default 64G)
#   ARF_BENCH_SLURM_TIME  time limit    (default 08:00:00)
#   ARF_BENCH_SBATCH_OPTS extra sbatch flags (e.g. "--partition=batch")
# Any ARF_BENCH_* sweep variable set in the environment still applies on top
# (e.g. ARF_BENCH_ITERS=3).
set -euo pipefail
cd "$(dirname "$0")/.."

CPUS="${ARF_BENCH_SLURM_CPUS:-17}"
# Nodes have ~1TB RAM; expct's cartesian hit the old 64G ceiling and OOM'd at
# n=20000. Request generously (jobs run on separate nodes / spread out anyway).
MEM="${ARF_BENCH_SLURM_MEM:-256G}"
TIME="${ARF_BENCH_SLURM_TIME:-08:00:00}"
OPTS="${ARF_BENCH_SBATCH_OPTS:-}"
mkdir -p bench/results bench/logs

# name|env-settings (applied per job on top of the caller's environment)
CONFIGS=(
  "forde|ARF_BENCH_OPS=forde"
  "lik|ARF_BENCH_OPS=lik"
  "adversarial-rf|ARF_BENCH_OPS=adversarial_rf"
  "forge-sep-small|ARF_BENCH_OPS=forge ARF_BENCH_ROWMODE=separate ARF_BENCH_NEVIDENCE=100 ARF_BENCH_NSYNTH=1"
  "forge-sep-large|ARF_BENCH_OPS=forge ARF_BENCH_ROWMODE=separate ARF_BENCH_NEVIDENCE=1000 ARF_BENCH_NSYNTH=100"
  "forge-or-small|ARF_BENCH_OPS=forge ARF_BENCH_ROWMODE=or ARF_BENCH_NEVIDENCE=100 ARF_BENCH_NSYNTH=1"
  "forge-or-large|ARF_BENCH_OPS=forge ARF_BENCH_ROWMODE=or ARF_BENCH_NEVIDENCE=1000 ARF_BENCH_NSYNTH=100"
  "expct-sep|ARF_BENCH_OPS=expct ARF_BENCH_ROWMODE=separate ARF_BENCH_NEVIDENCE=100"
  "expct-or|ARF_BENCH_OPS=expct ARF_BENCH_ROWMODE=or ARF_BENCH_NEVIDENCE=100"
)

for cfg in "${CONFIGS[@]}"; do
  name="${cfg%%|*}"
  env_settings="ARF_BENCH_LABEL=$name ${cfg#*|}"
  cmd=(sbatch --job-name="arf-bench-$name"
       --cpus-per-task="$CPUS" --mem="$MEM" --time="$TIME"
       --output="bench/logs/slurm-$name-%j.log")
  # shellcheck disable=SC2206  # intentional word-splitting of extra flags
  [ -n "$OPTS" ] && cmd+=($OPTS)
  cmd+=(--wrap="ARF_BENCH_SCRIPT=bench/sweep-ops.R $env_settings bench/run-sweep.sh clean")
  if [ "${DRY_RUN:-0}" = "1" ]; then
    printf '%q ' "${cmd[@]}"; echo
  else
    "${cmd[@]}"
  fi
done
