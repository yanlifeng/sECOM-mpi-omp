#!/usr/bin/env bash
# Run only in a new p3 scratch directory, under an approved taskset 24-31.
set -eu -o pipefail
test ! -e test.log && test ! -e check-runtime || exit 91
test "$(awk '/Cpus_allowed_list/{print $2}' /proc/self/status)" = 24-31 || exit 92
module load mpi/openmpi-x86_64
{
date -Ins
hostname
mpifort --showme
mpirun --version | head -1
grep Cpus_allowed_list /proc/self/status
sha256sum secom_mpi.f90 extrnl.f mpideck check_packed.f90 runtime_place.c runtime_rank.sh prepare.py halo9.f90.inc native_unit.sh
for cpu in {24..31}; do
    p=/sys/devices/system/cpu/cpu$cpu/topology
    echo "$(cat "$p/physical_package_id"):$(cat "$p/core_id") $(cat "$p/thread_siblings_list")"
done
flags=(-fopenmp -fcheck=all -O3 -g -fno-fast-math -ffp-contract=off -ffpe-trap=invalid,zero,overflow)
mpicc -O2 -c runtime_place.c
mpifort "${flags[@]}" -cpp -DSECOM_USE_MPI -c secom_mpi.f90
mpifort "${flags[@]}" check_packed.f90 secom_mpi.o runtime_place.o -o check-runtime
sha256sum check-runtime
for args in '5 1 1' '5 2 1' '7 4 1' '4 4 1' '7 2 2' '5 1 2' '11 4 2'; do
  set -- $args
  echo "BEGIN JM=$1 ranks=$2 OMP=$3"
  timeout 15s mpirun --host localhost:4 -np "$2" --nooversubscribe --bind-to none \
      --mca plm isolated --mca pml ob1 --mca btl self,vader sh runtime_rank.sh "$@"
  echo PASS
done
date -Ins
} 2>&1 | tee test.log
sha256sum test.log
