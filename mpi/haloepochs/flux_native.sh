#!/usr/bin/env bash
# From a fresh directory, within an explicitly allocated four-core subset.
# taskset -c <LO>-<HI> bash /overlay/flux_native.sh /frozen/parent/source
set -eo pipefail
test "$#" = 1
src=$(cd "$1" && pwd)
kit=$(cd "$(dirname "$0")" && pwd)
test ! -e unit && test ! -e test.log || exit 91
mask=$(awk '/Cpus_allowed_list/{print $2}' /proc/self/status)
[[ "$mask" =~ ^([0-9]+)-([0-9]+)$ ]] || exit 92
first=${BASH_REMATCH[1]}; last=${BASH_REMATCH[2]}
(( first >= 4 && last <= 35 && last-first == 3 )) || exit 92
module load mpi/openmpi-x86_64
{
  python3 -B "$kit/flux_direction.py" --source "$src" --out unit
  cd unit
  sha256sum secom_mpi.f90 advave.f extrnl.f "$kit/check_flux_dir.f90" "$src/mpideck"
  mpifort --version | head -1
  ff=(-cpp -DSECOM_USE_MPI -O3 -fopenmp -std=legacy -fcheck=all -ffixed-line-length-132
      -ffree-line-length-none -fno-fast-math -ffp-contract=off -ffpe-trap=invalid,zero,overflow
      -march=native -mprefer-vector-width=256 -I"$src")
  mpifort "${ff[@]}" -c secom_mpi.f90 extrnl.f advave.f
  mpifort "${ff[@]}" "$kit/check_flux_dir.f90" secom_mpi.o -o checkflux
  export OMP_PROC_BIND=TRUE OMP_DYNAMIC=FALSE
  for spec in '1 2 5 none' '2 1 7 none' '4 1 4 none' '4 1 11 none' '2 2 5 none' '2 2 7 reverse'; do
    set -- $spec
    rankc=$1; nt=$2; ny=$3; mutation=$4; log="case-${rankc}-${nt}-${ny}-${mutation}.log"
    echo "TEST ranks=$rankc threads=$nt JM=$ny mutation=$mutation mask=$mask"
    rc=0
    OMP_NUM_THREADS=$nt timeout 20s mpirun --host localhost:4 --nooversubscribe --bind-to none \
      --mca plm isolated --mca pml ob1 --mca btl self,vader -np "$rankc" \
      bash -c 'l=$(($1+OMPI_COMM_WORLD_LOCAL_RANK*OMP_NUM_THREADS)); h=$((l+OMP_NUM_THREADS-1))
        export OMP_PLACES="{$l}"; if [ "$h" != "$l" ]; then OMP_PLACES+=",$(printf "{%s}" "$h")"; fi
        exec taskset -c "$l-$h" ./checkflux "$2" "$3" "$4" "$5"' \
      rank "$first" "$rankc" "$nt" "$ny" "$mutation" > "$log" 2>&1 || rc=$?
    if [ "$mutation" = none ]; then
      test "$rc" = 0 && grep 'PASS_FLUX_DIR' "$log" || { cat "$log"; exit 90; }
    else
      test "$rc" = 41 && grep 'DETECTED_WRONG_DIRECTION' "$log" || { cat "$log"; exit 90; }
    fi
  done
} 2>&1 | tee test.log
