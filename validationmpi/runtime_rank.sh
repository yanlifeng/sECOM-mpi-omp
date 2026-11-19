#!/bin/sh
# Test fixture only; CPU24..31 lie inside the p3 CPU4..35 budget.
set -eu
base=$((24 + 2 * OMPI_COMM_WORLD_LOCAL_RANK))
case "$3" in
  1) cpus=$base; OMP_PLACES="{$base}";;
  2) cpus="$base,$((base+1))"; OMP_PLACES="{$base},{$((base+1))}";;
  *) exit 94;;
esac
export OMP_PLACES OMP_NUM_THREADS="$3" OMP_PROC_BIND=close OMP_DYNAMIC=false
exec taskset -c "$cpus" ./check-runtime "$@"
