# 澳门天鸽 sECOM OpenMP 研究版 v6

这是从附件和已验证的 v2 衍生的单机 CPU 版本，不是上游 MPI。包内 `source/` 是完整源码；`case20m/input/` 是 24 个冻结输入，保留全网格、4,800 个内模步、每步 10 个外模子步、四个平均场和最终重启。不要在输入模板里直接运行。

有 GNU Fortran、OpenMP 和 netCDF-Fortran 的机器，从解压目录执行：

```sh
make -C source -f Makefile.portable -j4 MODE=serial
make -C source -f Makefile.portable -j4 MODE=omp
python3 smoke.py /new/path/secom-smoke 4
```

若使用 p3 的 GCC 8.5，将下面两个变量作为 make 命令行参数传给两个构建：

```sh
NF_CONFIG=/home/user_home/ylf/kt5-licom-secom-ylf/deps/bin/nf-config
FFLAGS='-O3 -fno-fast-math -ffp-contract=off -std=legacy -ffixed-line-length-132 -g -fbacktrace -ffpe-trap=invalid,zero,overflow'
```

例如：`make -C source -f Makefile.portable -j4 MODE=omp NF_CONFIG="$NF_CONFIG" FFLAGS="$FFLAGS"`。p3 的 HDF5 库需要相同的初始化兼容层：串、并行都使用已记录 SHA 的 `h5_preinit.so`（包内附 C 源码）；设置相应的 `LD_LIBRARY_PATH`。这是库兼容项，不是加速手段，也没有放宽求解器浮点检查。

正式运行：将 `case20m/input` 复制到全新目录，进入副本，运行解压包 `source/build/omp/s--pom.x` 的绝对路径。在 p3 设置 `OMP_NUM_THREADS=32 OMP_DYNAMIC=FALSE OMP_PROC_BIND=TRUE`，并令 `OMP_PLACES` 为物理核 40 至 71 的 32 个单核集合；在其他机器先核对拓扑，不照搬 CPU 编号。

同包 `MODE=serial` 或 1 线程运行包含 v6 的算法整理，**不是历史原串行性能基准**。本次两种加速比均以进程启动到退出计时，包含实际读写。不要把结果外推到全台风、三维、MPI、其他机器或实测淹没精度。未改变输出频率和产品，不使用 fast-math。
