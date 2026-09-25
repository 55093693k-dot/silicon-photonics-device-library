# 仿真设置说明 — `PSR_ROTATE_ADIABATIC_BLT_v1`

> 配套工程文件：`run_bilevel_psr.py`（本文件夹内，可直接运行）｜ 求解器：**Tidy3D EME（3D）**，tidy3d 2.12
> 数据与结果：`../data/S_PARAMS.csv`、`../data/PORT_MODE_NEFF.csv`、`../data/PORT_MODE_FIELDS.png`

## 1. 运行方式

```powershell
# 0) 依赖：pip install tidy3d==2.12；配置 API key（~/.config/tidy3d/config.toml）
python run_bilevel_psr.py --stage geometry      # 出整器件结构图（本地，免费）
python run_bilevel_psr.py --stage modes         # 本地模式复核（免费）：TE0/TM0 有效折射率随 x 的演化
python run_bilevel_psr.py --stage eme           # 提交前校验 + 打印**估算成本**（不提交，免费）
python run_bilevel_psr.py --stage eme --submit  # 真正提交（按估计成本计费，见 §6）
```
脚本在 `--stage eme` 时会先打印**几何核对报告**并调用 `validate_pre_upload()`，通过后才提交。

## 2. 结构域与坐标

| 项 | 值 |
|---|---|
| 传播方向 | **+x**（`EMESimulation(axis=0)`） |
| 横向分离方向 | **y**（两条波导在 y 方向分开；仿真域 y = ±6 µm 以完整包住） |
| 高度 | z（Si 220 nm + 90 nm 部分刻蚀平板，单层硅） |
| 仿真域 | `center = (248.916, 0, 0.11) µm`，`size = (537.832, 12.0, 1.75) µm` |
| 器件 x 范围 | 0 … **537.832 µm**（= L_blt 100 + L_s 5 + L_ac 300 + L_t 30 + 2Rθ） |

## 3. 网格（Mesh）

| 项 | 设置 |
|---|---|
| 横截面网格 | `GridSpec.auto(min_steps_per_wvl=20, wavelength=1.55 µm)` |
| 传播方向网格（EME 单元划分） | **`EMEExplicitGrid`（157 个单元，156 个内部边界）**，在几何突变处对齐：<br>锥段 0→100 µm **2 µm/格**；突变界面区 95→110 µm **0.5 µm/格**；耦合器与输出 110→537.8 µm **5 µm/格** |
| 每单元模式规格 | `EMEModeSpec(num_modes=2)` → **端口模式数 = 2**（TE0、TM0） |
| 网格收敛性 | 单元数 60 → 157 时转换效率最小值 0.709 → 0.924 并趋于平坦（已收敛方向正确） |

## 4. 边界条件

- **EME 不使用 PML/周期边界**：两端以**端口**截断，端口相对仿真域边界内缩 `port_offsets = (10, 10) µm`。
- 因此端口平面位于 **x = 0 µm（输入）** 与 **x = 537.832 µm（输出）**，端点处均为**z 向开放**的介质界面。

## 5. 模式源 / 端口定义

| 项 | 设置 |
|---|---|
| 端口平面方向 | **x-normal**（垂直于传播方向 +x）；监视器 `size = (0, ∞, ∞)`，`center = (−10, 0, 0.11)`（输入）与 `(x_dev+10, 0, 0.11)`（输出） |
| 源定义 | **端口模式注入**（EME 端口模式即源）：输入端口模式 0 = **TE0**、模式 1 = **TM0**（按 neff 识别：2.412 / 1.809 @1.58 µm） |
| 输出端模式 | 0 = 宽臂（0.65 µm）**TE0**（neff 2.649）、1 = 窄臂（0.50 µm）**TE0**（2.498） |
| 端口模式求解 | `td.EMEModeSpec(num_modes=2)`；波导模式在端口截面处求解（含 90 nm 平板） |
| 约束 | `constraint="unitary"`（**注意：见 §7 已知事项**） |
| 端口模式存储 | `store_port_modes=True`（便于判定 mode_index ↔ 波导） |
| 附加监视器 | `EMEModeSolverMonitor`：`modes_in`、`modes_out`、`port_modes`；`EMECoefficientMonitor`：`coeffs` |

## 6. 波长与成本

| 项 | 值 |
|---|---|
| 频点 | **9 个波长：1.50 → 1.58 µm**（步长 0.01 µm） |
| 预计成本 | **≈1.335 FlexCredit**（2 模式 / 157 单元；4 模式 / 157 单元 ≈1.573） |
| 输出数据 | 云端任务的**大体积结果文件**（hdf5 / 本地工程文件，单文件 19.2–57.4 MB）落在**本机临时目录**，**未纳入本仓**；需要离线重分析时按 §7 重跑一次 |

## 7. 已知事项（复现时请留意）

1. `constraint="unitary"` 会强制 S 矩阵酉性；当端口模式数很少（如 2）时，功率几乎全部留在这些模式内
   → **列和恰为 1、串扰量级极小**。此时**只引用相对分布与转换效率**，不要引用绝对插损/回损。
2. 端口模式数取 4 时，高阶模式接近包层截止（neff≈1.44）并出现**频率间模式跳变**，会使某些 S 元素非物理（`|S21|²>1`）；
   取 2 可获得物理自洽的结果。
3. 材料模型：Si 与 SiO₂（材料库模型；介质定义见脚本内 `si_medium()` 与 `td.material_library["SiO2"]`）。
4. 长任务请用后台进程运行并轮询日志；stdout 需 `python -u`（脚本内含 UTF-8 输出处理）。

## 8. 数据依赖（自洽性）

| 项 | 位置 | 说明 |
|---|---|---|
| 仿真工程文件（**本仓自足**） | `sim/`：`run_bilevel_psr.py`、`run_eme_coupler.py`、`run_fdtd_short.py`、`plot_port_mode_fields.py` | **不依赖仓外的私有文件**：几何 / 材料 / 网格 / 源 / 端口定义全部写在脚本内 | 
| 验证数据（**本仓内**） | `../data/S_PARAMS.csv`、`PORT_MODE_NEFF.csv`、`PORT_MODE_FIELDS.png` | 可直接与自己跑出的结果逐点对照 | 
| 分析输入 hdf5（**大体积，未纳入本仓**） | 云端任务输出，单文件 19.2–57.4 MB（合计 ≈163 MB） | **是结果文件而不是脚本依赖**（可由 `sim/` 重跑生成，需云额度）⇒ 想 0 成本重分析，需要先自己重跑一次 |
| 本地求解器工程（**未纳入本仓**） | `.lms`（1.60 GB）、`.fsp`（185 MB ×2） | 几何 / 材料 / 网格定义全在 `sim/*.py` 内 ⇒ **可完整重建** |

> 维护提醒：本仓**不含 hdf5**；如需引用离线分析结果，请按上表自己重跑一次再对照。
