# jinr 探测器程序说明（单机箱 XIA）

本文件只说明本次移植中**为 jinr 探测器新增或修改**的程序及其用法。
VME 相关程序、sift 接口以及其余探测器程序均原样保留，未在本文件中展开。

## 1. 概述

本次移植将原 RIBLL2026 三机箱 XIA 管线改到 JINR 单机箱 XIA 采集系统上：

- 只读取一个 XIA 机箱（crate 0，module 0–8，sid 2–10）。
- 在 `crush` 阶段完成原始二进制 decode 与探测器 mapping。
- 在 `coke` 阶段以 main_trigger 为时间轴锻造 `jinr_trigger` 与 `jinr_t1_trigger`。
- 在 `smelt` 阶段重构 `jinr_t0d1`、`jinr_t0d2`、`jinr_beam`、`jinr_gagg` 事件。

各 jinr 探测器与输出文件的对应关系：

| 探测器 | 含义 | 生成阶段 | 输出文件 |
|--------|------|----------|----------|
| jinr_trigger | main trigger（module 0 ch2） | crush → coke | `grit/trigger_%04d.root` → `ingot/trigger_%04d.root` |
| jinr_t1_trigger | t1 trigger（slot 内出现 module 0 ch3 vme_trigger） | coke | `ingot/trigger_t1_%04d.root` |
| jinr_t0d1 | T0D1 DSSD（module 0/1/2） | crush → smelt | `grit/t0d1_%04d.root` → `ingot/t0d1_%04d.root` |
| jinr_t0d2 | T0D2 DSSD（module 3/4/5/6/8） | crush → smelt | `grit/t0d2_%04d.root` → `ingot/t0d2_%04d.root` |
| jinr_beam | TOF（module 0 ch15） | crush → smelt | `grit/beam_%04d.root` → `ingot/beam_%04d.root` |
| jinr_gagg | GAGG（module 7） | crush → smelt | `grit/gagg_%04d.root` → `ingot/gagg_%04d.root` |

> 说明：探测器概念名使用 `jinr_` 前缀，实际 ROOT 文件名沿用原有短名（`t0d1`、`t0d2`、`beam`、`gagg`、`trigger`），未额外加前缀。

## 2. 编译方法

依赖：

- CMake 3.10 或更高
- 支持 C++17 的编译器
- CERN ROOT 6.20 或更高（通过 `find_package(ROOT ... CONFIG REQUIRED)`）

本项目要求 out-of-source 编译。在仓库根目录执行：

```bash
cmake -S . -B build
cmake --build build -j
```

编译完成后，本文件涉及的 jinr 相关可执行文件位于：

| 程序 | 可执行文件路径 |
|------|----------------|
| crush_trigger / crush_beam / crush_t0d1 / crush_t0d2 / crush_gagg | `build/bin/crush/` |
| coke | `build/bin/coke` |
| smelt | `build/bin/smelt` |

即：

```
build/bin/crush/crush_trigger
build/bin/crush/crush_beam
build/bin/crush/crush_t0d1
build/bin/crush/crush_t0d2
build/bin/crush/crush_gagg
build/bin/coke
build/bin/smelt
```

> 若 ROOT 未安装在默认路径，可先执行 `source <ROOT安装路径>/bin/thisroot.sh` 设置环境后再运行 cmake。

## 3. 文件路径约定

- 原始 XIA 数据（crate 0）：
  `<paths.raw_xia0>/<run:%04d>/<prefix.raw_xia0>_R<run:%04d>_M<module:%02d>.bin`
- crush 输出（grit）：
  `<workspace>/<paths.grit>/<detector>_<run:%04d>.root`
- coke / smelt 输出（ingot）：
  `<workspace>/<paths.ingot>/<detector>_<run:%04d>.root`（main）
  `<workspace>/<paths.ingot>/<detector>_t1_<run:%04d>.root`（t1）

以上目录与文件名前缀由 `config.toml` 控制，默认 `-c config.toml`。

## 4. crush 阶段（decode + mapping）

每个探测器一个独立程序，均以 run 号作为位置参数，可选 `-c` 指定配置文件：

```
crush_trigger <run> [-c config.toml]
crush_beam    <run> [-c config.toml]
crush_t0d1    <run> [-c config.toml]
crush_t0d2    <run> [-c config.toml]
crush_gagg    <run> [-c config.toml]
```

各程序读取的 module 与采样率：

| 程序 | module | 采样率（MHz） | 映射要点 |
|------|--------|---------------|----------|
| crush_trigger | 0 | 500 | ch2=Main trigger，ch3=Vme trigger，带 external_time |
| crush_beam | 0 | 500 | ch15=TOF |
| crush_t0d1 | 0, 1, 2 | 500, 500, 100 | m0-ch0=b15；m1=f0–15；m2=b0–14 |
| crush_t0d2 | 3, 4, 5, 6, 8 | 100 | 含 f22/b1/b23/b31 归属与 error 通道修补 |
| crush_gagg | 7 | 100 | index=15−ch，跳过 ch6（apd10 不存在），编号 0–13 |

## 5. coke 阶段（trigger 锻造）

```
coke --run <run> [-c config.toml]
# 或
coke --xia-run <run> [-c config.toml]
```

- 输入：`grit/trigger_%04d.root`
- 输出：
  - `ingot/trigger_%04d.root`：jinr_trigger（main trigger）
  - `ingot/trigger_t1_%04d.root`：jinr_t1_trigger（slot 内出现 ch3 vme_trigger 的 main trigger）

以 main_trigger 为时间轴，采用环形缓冲 slot + 最近邻匹配（window=1000 ns，kSlotNum=1000）聚合事件。

## 6. smelt 阶段（探测器事件重构）

```
smelt --run <run> [-t main|t1] <detector...>
```

`<detector...>` 为位置参数，可一次指定多个。本次 jinr 相关可用值：`t0d1`、`t0d2`、`beam`、`gagg`。

- 输入：`ingot/trigger_%04d.root`（`-t main`）或 `ingot/trigger_t1_%04d.root`（`-t t1`）
- 输出：
  - `-t main`：`ingot/<detector>_%04d.root`
  - `-t t1`：`ingot/<detector>_t1_%04d.root`

> 当前 jinr 探测器暂只提供 main 版本（`-t main`）；t1 版本待 sift 引入 VME 事件后再给出。

示例（锻造 main 版本）：

```
smelt --run 123 t0d1 t0d2 beam gagg
```

## 7. 完整运行示例

以 run 123 为例，按顺序执行：

```bash
# 1) decode + mapping（crush）
crush_trigger 123
crush_beam    123
crush_t0d1    123
crush_t0d2    123
crush_gagg    123

# 2) trigger 锻造（coke）
coke --run 123

# 3) 探测器事件重构（smelt，main 版本）
smelt --run 123 t0d1 t0d2 beam gagg
```