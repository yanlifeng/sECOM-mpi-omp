#!/usr/bin/env bash
# From a fresh directory on an explicitly assigned contiguous four-core slot:
# taskset -c LO-HI bash native_health.sh R007/source check_extrema_mpi.f /kit/extrema2/make_fixture.py
set -eo pipefail
ulimit -c 0
test "$#" = 3
kit=$(cd "$(dirname "$0")" && pwd)
src=$(cd "$1" && pwd)
orig=$(realpath "$2"); mkfixture=$(realpath "$3")
test ! -e qa && test ! -e test.log || exit 91
mask=$(awk '/Cpus_allowed_list/{print $2}' /proc/self/status)
[[ "$mask" =~ ^([0-9]+)-([0-9]+)$ ]] || exit 92
first=${BASH_REMATCH[1]}; last=${BASH_REMATCH[2]}
(( first >= 4 && last <= 35 && last-first == 3 )) || exit 92
module load mpi/openmpi-x86_64
{
  python3 -B "$kit/prepare.py" --source "$src" --out qa
  python3 -B "$mkfixture" "$orig" qa/expanded.f
  cd qa
  sha256sum secom_mpi.f90 baro_check.f expanded.f "$kit/check_health.f" "$orig" "$src/mpideck"
  mpifort --version | head -1
  ff=(-cpp -DSECOM_USE_MPI -O3 -fopenmp -std=legacy -fcheck=all -ffixed-line-length-132
      -ffree-line-length-none -fno-fast-math -ffp-contract=off -ffpe-trap=invalid,zero,overflow
      -g -fbacktrace -march=native -mprefer-vector-width=256 -I"$src")
  mpifort "${ff[@]}" -c secom_mpi.f90 baro_check.f
  mpifort "${ff[@]}" "$orig" secom_mpi.o baro_check.o -o stock
  mpifort "${ff[@]}" expanded.f secom_mpi.o baro_check.o -o extended
  mpifort "${ff[@]}" "$kit/check_health.f" secom_mpi.o baro_check.o -o contract
  export OMP_PROC_BIND=TRUE OMP_DYNAMIC=FALSE
  for exe in stock extended contract; do
    for conf in '2 1' '4 1' '2 2'; do
      set -- $conf; nr=$1; nt=$2
      echo "TEST $exe ranks=$nr nt=$nt mask=$mask"
      second=$nt; if [[ "$exe" = contract ]]; then second=ok; fi
      OMP_NUM_THREADS=$nt timeout 20s mpirun --host localhost:4 -np "$nr" --nooversubscribe --bind-to none \
        --mca plm isolated --mca pml ob1 --mca btl self,vader bash -c \
        'l=$(($1+OMPI_COMM_WORLD_LOCAL_RANK*OMP_NUM_THREADS)); h=$((l+OMP_NUM_THREADS-1)); export OMP_PLACES="{$l}"
        if [ "$h" != "$l" ]; then OMP_PLACES+=",$(printf "{%s}" "$h")"; fi
        exec taskset -c "$l-$h" "./$2" "$3" "$4"' rank "$first" "$exe" "$nr" "$second"
    done
  done
  for mode in trap-old trap-new; do
    rc=0
    OMP_NUM_THREADS=1 timeout 20s mpirun --host localhost:2 -np 2 --nooversubscribe --bind-to none \
      --mca plm isolated --mca pml ob1 --mca btl self,vader bash -c \
      'l=$(($1+OMPI_COMM_WORLD_LOCAL_RANK)); export OMP_PLACES="{$l}"; exec taskset -c "$l" ./contract 2 "$2"' \
      rank "$first" "$mode" > "$mode.log" 2>&1 || rc=$?
    [[ $rc != 0 && $rc != 124 ]] && grep -q EXPECT_WINNER_SIGFPE "$mode.log" &&
      grep -Eqi 'received signal SIGFPE|floating.point exception|signal 8 ' "$mode.log" &&
      ! grep -q 'STOP 87' "$mode.log" || { cat "$mode.log"; exit 87; }
    echo "EXPECTED selected speed row SIGFPE: $mode rc=$rc"
  done
} 2>&1 | tee test.log
