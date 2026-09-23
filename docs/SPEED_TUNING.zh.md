# 在新集群上为 GOCIA 调优 VASP / CP2K 弛豫速度

英文版 `SPEED_TUNING.md`。工具在 `tools/speedtest/`。参考结果见第 7 节。

GOCIA 几乎所有时间都花在子代结构的局部弛豫上。要最大化的量是**给定核数预算下每小时弛豫完的结构数, 且精度不变**。一个设置只有在更快**并且**通过第 4 节精度门时才采用。下面的方法与调度系统、集群无关; 唯一和集群相关的文件是 `tools/speedtest/` 的 JSON 配置。

## 1. 需要准备

- 生产流程的**分级模板**: VASP 的 `INCAR-1..3`, 或 CP2K 的 `cp2k-1..3.inp` (GOCIA 模板约定, 用 `@INCLUDE` / `@SET`)。
- 生产体系的**一个已弛豫结构** (每步测试, 第 3 节 B-G)。
- 同一体系**至少 4 个未弛豫的 GA 子代**: 2 个变异子代 (如 rattle) 和 2 个交叉子代, 即真实的远离平衡的起点 (每结构测试, H-I)。每次 GA 都会碰到需要长距离重排的结构; 手上有就留一个, 优化器对比由它决定。
- 调度系统信息: 节点类型 (CPU 型号、每节点核数与内存)、每用户核数上限、小作业是否共享节点、整节点作业是否排很久。

## 2. 测量规则

1. **指标 = 墙钟秒/力计算** (作业墙钟 ÷ 力计算次数) 用于每步测试; **墙钟分钟/结构** (各级之和) 用于弛豫。程序内部计时只用于诊断; 若要用, 先在一个运行上与墙钟对账 (差 >5% 说明计时漏了一段)。两种代码的力计算次数口径一致: VASP `NSW=5` 为 5 次, CP2K `MAX_ITER 5` 为 6 次。
2. **变体放在同一次分配里比较。** 共享节点上邻居负载让运行时间波动 10–30%。`speedtest.py ab` 把所有变体放进同一个作业、轮换顺序运行; 至少 2 个作业 × 2 轮 (每个变体 ≥4 个样本), 看**作业内相对参照变体的比值**。
3. **用确定性指标区分原因与噪声**: 每次力计算的 SCF 迭代数、首步能量 (同一几何)、N 步后的终态能量与几何, 都不受节点负载影响。
4. **计时用单独作业, 不用并行打包的运行。** 若多个运行在同一分配里并行, 要逐个核对绑核 (`ps -o psr`); 这里出错会让所有运行挤在同一批核上, 看起来像负载高。
5. **每个运行记录所在主机**, 保留原始输出; 没有主机和输入的结果无法复用。
6. **按共址数值做规划。** GOCIA 是很多 worker 并排跑; 同时提交 10–12 个相同作业, 用它们的均值, 不用单独运行的最好成绩。

## 3. 步骤

每一步: 变什么、测什么、何时采用。B–G 用已弛豫结构做 5 步几何优化; H–I 用子代。

**A. 盘点与编译。** 列出节点类型与限制; 编译或找到每个程序。CP2K: 针对目标 CPU 编译 (对应微架构的编译选项), 能带 ELPA 就带; 用回归测试确认。计时之前先对照第 5 节的已知编译/运行陷阱。

**B. 节点类型 × 核数。** 每种节点上单独作业 8、16、24、48 (96) 核, 两种代码, 纯 MPI。选最快的节点类型, 再选相对 8 核并行效率仍 ≥90% 的最大核数; 超过它, 每结构多用的核会降低吞吐。Γ 点板模型用很大的 rank 数还可能 SCF 不收敛; 不要超出实测范围。

**C. 共址。** 同时跑 10–12 个选定规模的相同作业; 其均值与单独运行之比。所有规划都用共址数值。

**D. 并行布局。**
- VASP: `NCORE` 取 rank 数的各个因子 (√ranks 只是起点), `NSIM` 4/8/16, 只有多个 k 点时才用 `KPAR`, MPI × OpenMP (如 8×2、4×4)。
- CP2K: MPI × OpenMP (1/2/4 线程)。

**E. 库与编译版本。** 同样输入换程序或库: 对角化库 (ScaLAPACK 与 ELPA, `&GLOBAL PREFERRED_DIAG_LIBRARY`)、通用编译与针对 CPU 的编译、AMD 上带 AVX-512 时 MKL 的代码路径 (`tools/speedtest/mkl_amd_shim.c`)。能量必须完全一致 (≤0.01 meV) 才采用。

**F. SCF 设置, 分两种情形测。**
- CP2K: `EPS_SCF` (1e-5 / 1e-4 / 1e-3)、Broyden 混合里的 Kerker 阻尼 `BETA` (1.0 / 1.5 / 2.5)、`ALPHA`、`NBROYDEN`、Pulay 混合、波函数外推、`ADDED_MOS`、FFTW 规划方式。金属用对角化 + Fermi-Dirac smearing (OT 在金属板上收敛到与极小化器有关的整数占据态)。
- VASP: `EDIFF`、`ALGO` (Fast / VeryFast / Normal), 流程允许时才动 `NELMIN`。
- 情形 1: 已弛豫结构, 5 步 (靠近极小的小步长)。
- 情形 2: 一个子代, 第一级设置下的前 10 步 (大步长)。每步几何变化大时, 放宽 SCF 省得少得多, 有些设置表现也不同。同一个设置可能适合第一级而不适合第二、三级。

**G. 截断能与网格。** 三个结构在生产截断能与一个高截断能参照下做单点, 比较**相对**能量 (绝对能量不是判据)。改动 FFT 网格选项也要做同样检查 (第 5 节第 7 条)。

**H. 各级优化器。** 用全部子代跑完整的多级弛豫 (`speedtest.py pilot`), 每次只改一级的优化器: CP2K BFGS (信任半径)、LBFGS、CG; VASP CG (`IBRION 2`)、准牛顿 (`IBRION 1`)、VTST FIRE / LBFGS (`IBRION 3`, `IOPT 7` / `1`)。记录每级步数与墙钟, 并检查第一级能量是否下降 (接受上坡步的优化器会把预弛豫变成随机扰动)。每结构耗时通常由优化器决定, 而不是由任何每步设置决定。

**I. 最终试点与吞吐。** 全部选定设置一起跑全部子代; 每结构平均分钟数; 吞吐 = floor(核数上限 ÷ 每作业核数) × 60 ÷ 平均分钟 (`speedtest.py throughput`)。两个代码在同一组结构上比较, 并先确认两者到达了同一类极小 (第 4 节最后一条), 才能说一个代码更快。

## 4. 精度门

| 改动 | 判据 |
|---|---|
| 编译 / 库 / MKL 路径 | 首步与终态能量完全一致 (≤0.01 meV) |
| SCF 容差、混合 (第二、三级) | 同一几何首步 ΔE < 1 meV; 从极小出发 5 步后 ΔE_final < 1 meV 且最大位移 < 1e-3 Å |
| SCF 容差 (第一级预弛豫) | 同一几何首步 ΔE < 10 meV |
| 截断能 / FFT 网格 | 相对能量与高截断能参照相差 5–10 meV 以内 |
| smearing (CP2K) | 没有 "Add more MOs for proper smearing" 警告 |
| 优化器, 完整弛豫 | 各级都收敛; 终态在同一个势能面上比较 (同一程序、网格、紧 SCF 单点); 同一子代在很小的数值改动下可落进相差约 1 eV 的不同极小, 所以不要用少数几个结构的终态能量给设置排序 |
| 代码之间 | 终态吸附质位点分配相同 (不同泛函的能量不可比) |

## 5. 实际踩过的坑

1. **计时范围不一致。** 解析脚本对一个代码累加 SCF 迭代时间、对另一个取整离子步时间, 前者被低估 8–18%。两边都用墙钟/力计算。
2. **打包运行重叠。** 在子 shell 里计算的核分片不前进, 所有运行都落在 0..N-1 号核上。用 `ps -o psr` 检查。
3. **Intel MPI 2021 的编译器包装会把环境变量 `FCFLAGS`/`CFLAGS` 加到自己的参数前面。** configure 会 export `FCFLAGS`, 于是 autoconf 的模块测试里 `-J` 出现两次 ("gfortran: Only one -J option allowed")。解决: 写小包装脚本, 先 `unset FCFLAGS FFLAGS CFLAGS CXXFLAGS` 再 exec 真正的 `mpif90` / `mpicc`。
4. **Intel MPI 2021 只带到 gfortran 11 的 `mpi.mod`**; gfortran ≥12 时包装器退回 Intel 编译的模块, `use mpi` 编译失败。解决: 加 `-I<impi>/include/gfortran/11.1.0` (这组模块 gfortran 13 能读)。
5. **预编译的 CP2K (conda) 可能要求比集群更新的 glibc**; 从源码编译。
6. **CP2K 2026.2 对角化重分布。** `cp_fm_redistribute_init` 与 `cp_fm_redistribute_work_finalize` 在赋值之后又用类型构造器重置重分布设置, 于是每次对角化只在一部分 rank 上跑, `&FM_DIAG_SETTINGS` 不起作用。用 `&GLOBAL &FM_DIAG_SETTINGS PRINT_FM_REDISTRIBUTE` 检查 (所有 rank 都应参与); 把重置挪到赋值之前即可修复 (实测 16 rank 对角化 235 → 114 s)。
7. **`EXTENDED_FFT_LENGTHS` 会改变网格。** 一组输入在 `&GLOBAL` 里有这个关键字、另一组没有, 绝对能量差 0.75 eV (250 Ry) 和 1.25 eV (350 Ry), 相对能量偏 30–60 meV。两个"相同"设置的运行首步就对不上时, 先逐行比较完整输入, 再去解释。
8. **CP2K BFGS 用作预弛豫。** 默认 0.25 Å 信任半径下, 它在 GA 子代上接受上坡步 (25 步内能量摆动约 10 eV)。LBFGS 17–24 次力计算收敛。靠近极小的各级用带模型 Hessian 的 BFGS 最快; LBFGS 在那里多花线搜索的力计算。
9. **VASP `IBRION 1` 离极小远时**会发散 (能量跳几百 eV), 接在 FIRE 级之后又要爬 86–342 步。三级都用 FIRE 是最快且稳定的选择。
10. **大步长时放宽 SCF 效果变小**: 第一级设置下每次力计算的 SCF 迭代从约 16 降到约 9, 但要把容差放到 1e-3; 只改混合约省 10%。
11. **共享集群上的整节点作业**可能排队数天, 而 16 核作业立刻开跑; `examples/gcga_cp2k/per-job/` 的单 worker 启动方式就是为此准备的。

## 6. 工具 (`tools/speedtest/`)

```bash
cp tools/speedtest/cluster.example.json cluster.json          # 调度头、模块、MPI 启动命令、程序路径
# 每步 A/B (B-G): 由分级模板和已弛豫结构生成 5 步基准输入
python tools/speedtest/speedtest.py bench-input --code cp2k --template examples/gcga_cp2k/cp2k-2.inp \
       --structure relaxed.vasp --out inputs/
python tools/speedtest/speedtest.py ab cluster.json --code cp2k --inputs inputs/ \
       --variants tools/speedtest/examples/variants_cp2k.json --out runs/scf --jobs 2 --rounds 2
for j in runs/scf/j*/job.sh; do sbatch $j; done
python tools/speedtest/speedtest.py parse-ab runs/scf --ref base --csv scf.csv
# 完整弛豫 (H-I): 原样模板 对比 第一级换回 BFGS
python tools/speedtest/speedtest.py pilot cluster.json --code cp2k --stages examples/gcga_cp2k \
       --structures kid1.vasp kid2.vasp kid3.vasp kid4.vasp --out runs/pilot_ref
python tools/speedtest/speedtest.py pilot cluster.json --code cp2k --stages examples/gcga_cp2k \
       --structures kid1.vasp kid2.vasp kid3.vasp kid4.vasp --out runs/pilot_bfgs1 \
       --edits tools/speedtest/examples/edits_stage1_bfgs.json
for j in runs/pilot_*/*/job.sh; do sbatch $j; done
python tools/speedtest/speedtest.py parse-pilot runs/pilot_* --csv pilots.csv
python tools/speedtest/speedtest.py throughput --minutes 24.3 --cap 1300 --cores-per-job 16
```

变体 = 对主输入的精确文本替换 (`edits`, 每个目标必须恰好出现一次)、追加行 (`append`), 以及可选的 `ranks`、`omp`、`exe`、`env` (换程序、`LD_PRELOAD` 等)。

## 7. 参考结果

硬件: 192 核 AMD EPYC 9655 节点 (Zen 5, AVX-512), 与其他作业共享, 每结构 16 核, 纯 MPI。体系: Cu(100) 6×6×4 板 + 4 CO + 10 H (162 原子, 底两层固定, Γ 点, 真空); VASP RPBE, CP2K PBE (GTH / MOLOPT-SR)。三级力判据 0.5 / 0.2 / 0.1 eV/Å。

**CP2K** (模板 `examples/gcga_cp2k/cp2k-1..3.inp`)

| 项 | 值 |
|---|---|
| 编译 | 针对 CPU 编译并带 ELPA; `&GLOBAL PREFERRED_DIAG_LIBRARY ELPA`; AMD 上 MKL 走 AVX-512 路径 |
| SCF | 对角化 + Fermi-Dirac 300 K, Broyden ALPHA 0.4 / NBROYDEN 8 + Kerker BETA 1.5, ADDED_MOS = max(30, ⌈N_atoms/2⌉) |
| 第一级 | LBFGS, CUTOFF 250 / REL_CUTOFF 40, EPS_SCF 1e-3, MAX_ITER 25, 写出波函数 |
| 第二、三级 | BFGS, 350 / 50 Ry, EPS_SCF 1e-4, 从上一级波函数续算 |
| 关闭 | `EXTENDED_FFT_LENGTHS`、OT、OpenMP |

**VASP** (`examples/gcga_vasp_tuned/INCAR-1..3`)

| 项 | 值 |
|---|---|
| 并行 | 16 个 MPI rank, `NCORE 8`, `NSIM 16`, AMD 上 MKL 走 AVX-512 路径; 不用 OpenMP |
| SCF | `ALGO Fast`, 各级 `EDIFF 1e-4` |
| 优化器 | 三级都用 FIRE (`IBRION 3`, `IOPT 7`, `POTIM 0`, `MAXMOVE 0.2`, `TIMESTEP 0.1`; 需要 VTST 版) |
| 续算 | 第二级写 WAVECAR/CHGCAR, 第三级读取 (`ISTART 1`, `ICHARG 1`) |

**实测** (每结构, 4 个子代的均值)

| | 调优前 | 调优后 | 1300 核每小时结构数 |
|---|---|---|---|
| CP2K | 95.8 min | 24.3 min | 200 |
| VASP | 155.6 min | 77.8 min | 62 |

每步手段 (5 步基准, 作业内比值): CP2K ELPA 0.90×, ELPA + MKL AVX-512 路径 0.80×, EPS_SCF 1e-4 0.64× (ΔE 0.06 meV), 再加 Kerker 1.5 为 0.57×; 单独针对 CPU 编译、减少 ADDED_MOS、FFTW MEASURE、REL_CUTOFF 40、DGEMM 网格后端、Pulay 混合和其他外推都没有收益或更慢。VASP NSIM 16 0.87×, MKL AVX-512 路径 0.90×, EDIFF 1e-4 0.90× (ΔE 0.41 meV); OpenMP 慢 1.7–2.3×, ALGO VeryFast 无收益。每结构手段: CP2K 第一级 LBFGS 0.67×; VASP 第二级 FIRE 的 SCF 工作量是 CG 的一半。最难的子代 (约 1 eV 的长下坡重排) CP2K 36 min、VASP 178 min, 两者到达同一吸附质构型。
