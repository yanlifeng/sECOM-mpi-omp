#!/usr/bin/env bash
# From a new scratch directory, within an explicitly assigned 4-core subset.
# Usage: taskset -c 20-23 bash native_unit.sh /snapshot/source /original/check_extrema_mpi.f
# First allowed CPU is discovered from /proc; no hardcoded claims on a core.
set -eo pipefail
test "$#" = 2
overlay=$(cd "$(dirname "$0")" && pwd)
parent=$(cd "$1" && pwd)
fixture=$(realpath "$2")
test ! -e unit && test ! -e test.log || exit 91
allowed=$(awk '/Cpus_allowed_list/{print $2}' /proc/self/status)
[[ "$allowed" =~ ^([0-9]+)-([0-9]+)$ ]] || exit 92
lo=${BASH_REMATCH[1]}; hi=${BASH_REMATCH[2]}
(( lo >= 4 && hi <= 35 && hi-lo == 3 )) || exit 92
module load mpi/openmpi-x86_64
{
  python3 -B "$overlay/prepare.py" --source "$parent" --out unit
  python3 -B "$overlay/make_fixture.py" "$fixture" unit/extended.f
  cd unit
  sha256sum baro_check.f secom_mpi.f90 extended.f "$fixture" "$parent/comdeck" "$parent/mpideck"
  mpifort --version | head -1
  export OMP_DYNAMIC=FALSE OMP_PROC_BIND=TRUE
  flags=(-cpp -DSECOM_USE_MPI -O3 -fopenmp -std=legacy -fcheck=all -g -fbacktrace
         -ffixed-line-length-132 -ffree-line-length-none -fno-fast-math -ffp-contract=off
         -ffpe-trap=invalid,zero,overflow -march=native -mprefer-vector-width=256 -I"$parent")
  mpifort "${flags[@]}" -c secom_mpi.f90 baro_check.f
  mpifort "${flags[@]}" "$fixture" secom_mpi.o baro_check.o -o stock
  mpifort "${flags[@]}" extended.f secom_mpi.o baro_check.o -o expanded
  for exe in stock expanded; do
    for conf in 2:1 4:1 2:2; do
      nr=${conf%:*}; nt=${conf#*:}
      echo "CHECK $exe ranks=$nr threads=$nt allowed=$allowed"
      OMP_NUM_THREADS=$nt timeout 20s mpirun --host localhost:4 --nooversubscribe --bind-to none \
        --mca plm isolated --mca pml ob1 --mca btl self,vader -np "$nr" \
        bash -c 'l=$(($1+OMPI_COMM_WORLD_LOCAL_RANK*OMP_NUM_THREADS)); h=$((l+OMP_NUM_THREADS-1))
          export OMP_PLACES="{$l}"; if [ "$h" != "$l" ]; then OMP_PLACES+=",$(printf "{%s}" "$h")"; fi
          exec taskset -c "$l-$h" "./$2" "$3" "$OMP_NUM_THREADS"' rank "$lo" "$exe" "$nr"
    done
  done
} 2>&1 | tee test.log
