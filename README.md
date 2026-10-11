# Timing-Driven MacroPlacer

ASPDAC 2026 时序驱动宏布局算法的源码仓库，基于 OpenROAD `mpl2`，包含全局宏布局和增量宏布局实现。

本仓库仅上传算法源码、源码补丁、构建配置和运行辅助脚本。论文 PDF、预编译二进制、基准数据、实验输出、日志、指标表及测试结果样本均未上传。

## 文件结构

- `01_source_code/mpl2/`：C++ 算法源码、头文件、Tcl/SWIG 接口、CMake 配置、原有文档和测试代码。
- `01_source_code/openroad_mpl2.patch`：相对于指定 OpenROAD 版本的完整算法补丁。
- `01_source_code/openroad_base_info.txt`：原始 OpenROAD 版本及子模块信息。
- `02_docker/`：Docker 构建文件及运行辅助脚本。
- `03_benchmarks/<design>/run_*.tcl`：设计对应的运行配置脚本；不含输入数据。
- `scripts/`：批量运行和指标采集程序；不含采集结果。

主要算法文件位于 `01_source_code/mpl2/src/`：

- `hier_rtlmp.cpp` / `hier_rtlmp.h`：全局宏布局。
- `IncreMacroPlacement.cpp`：增量宏布局。
- `MacroPlacementSolver.cpp` / `MacroPlacementSolver.h`：宏布局求解器。

Tcl 入口为 `rtl_macro_placer` 和 `incremental_macro_placer`。

## 构建

基础 OpenROAD 版本为 `v2.0-17941-g7fb347f37`，提交：

```text
7fb347f37d0068b4a9344d7e6ae1b668360e174e
```

`mpl2` 是 OpenROAD 的模块，需在完整 OpenROAD 源码和依赖环境中构建。Dockerfile 会获取该版本及其子模块，应用源码补丁并构建 OpenROAD。在装有 Docker 的 Linux 主机上，从仓库根目录执行：

```bash
bash 02_docker/build_docker.sh
```

构建后的程序位于容器内 `/opt/OpenROAD/build/src/openroad`。依赖由 `openroad/ubuntu22.04-dev:latest` 基础镜像提供，构建需要网络。基础镜像未锁定摘要，实际依赖版本可能随镜像更新变化。

## 运行

运行脚本保留了原始归档的目录约定。运行前需自行在 `03_benchmarks/` 下准备所需 LEF、Liberty、SDC、ODB 等输入，或修改 Tcl 脚本使用自己的数据路径；本仓库不提供这些数据。原始基准来源为 ChiPBench-D：<https://huggingface.co/datasets/MIRA-Lab/ChiPBench-D>。

输入准备好后，可运行：

```bash
bash 02_docker/run_global.sh bp_be
bash 02_docker/run_global_incre.sh bp_be
```

直接运行自己的 Tcl 流程时：

```bash
docker run --rm -v "$PWD":/workspace/archive aspdac2026-mpl2:archive \
  -no_init /workspace/archive/your_flow.tcl
```

批量脚本可通过 `ARCHIVE` 和 `OPENROAD` 指定本地路径；默认路径沿用开发环境。`run_global_result_incre.tcl` 需要先生成对应的全局布局输出。

运行产生的 `04_results/`、日志和指标文件已加入 `.gitignore`。测试目录保留测试程序，不包含测试数据和预期结果；若运行原有集成测试，需从对应 OpenROAD 版本补齐测试资源。

## 源码许可

保留原始源码中的 OpenROAD 及其他版权声明、许可条款。使用和再分发时请遵循相应文件中的许可说明。
