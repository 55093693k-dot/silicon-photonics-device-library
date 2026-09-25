# 仿真设置说明 — `PSR_SPLIT_TE_TM_DC_v1`（Tidy3D 3D FDTD）

> 配套工程文件（本文件夹内，可直接运行）：`run_tidy3d_splitter_field.py`（**当前主用**）、
> `run_tidy3d_splitter.py`、`run_tidy3d_splitter_sweep.py`、`run_tidy3d_ref_guide.py`
> 求解器：**Tidy3D 3D FDTD**（云端计费）｜ 数据与结果：`../data/`、`../PSR_SPLIT_TE_TM_DC_v1_report.md`

## 1. 运行方式

```powershell
# 依赖：pip install tidy3d；配置 API key（~/.config/tidy3d/config.toml）
python run_tidy3d_ref_guide.py                       # 参考直波导对照（模式序/链路自检）
python run_tidy3d_splitter_field.py --source TE      # 主用：TE0 注入，场积分法
python run_tidy3d_splitter_field.py --source TM      # TM0 注入（尚未完成）
python run_tidy3d_splitter_field.py --source BOTH    # 两偏振一次跑
python run_tidy3d_splitter_field.py --hdf5 ..\data\SPLIT_FIELD_TE.hdf5 --tag TE   # 仅离线重分析（0 成本；输入已随仓提供，见 §8）
```

## 2. 结构域与坐标

| 项 | 值 |
|---|---|
| 传播方向 | **+x**（x = -2 → 33.5 µm，含 stub） |
| 两臂分离方向 | **z**（两臂中心 z = ±0.30 µm；间隙 0.15 µm） |
| 高度方向 | **y**（Si 220 nm，y = 0 → 0.22 µm） |
| 仿真域 | `center = (cx, 0, 0)`，`size = (35.5, 2.8, 6.0) µm`（y、z 留足包层） |
| 接入 | wg2 侧 3 µm S-bend（横向偏移 1.35 µm）+ 两端 2 µm stub |

## 3. 网格

| 项 | 设置 |
|---|---|
| 全域网格 | `GridSpec(wavelength = 1.55 µm)` + **MeshOverride dl = (0.02, 0.01, 0.01) µm**（x 20 nm、y 10 nm、z 10 nm），覆盖器件区 y = ±0.25 µm、z = ±1.0 µm |
| 依据 | 本项目的网格收敛结论：截面 / 器件区加密到 **≤20 nm**（基准 **10 nm**）即已收敛（见 `../README.md` §4） |
| 收敛性 | 网格筛选以"TE 直通 >95%"为判据扫描；**未做过 `certify` 式收敛复跑** ⇒ 见 §7 已知事项 |

## 4. 边界条件

四面 **PML**（`BoundarySpec.all_sides(PML())`）。

## 5. 源与监视器

| 项 | 设置 |
|---|---|
| 源 | `ModeSource` @ x = -1.7 µm，入 wg1（z = -0.30），`size = (0, 1.0, 1.0) µm`，`direction = "+"`；<br>`GaussianPulse(freq0 = c/1.55 µm, fwidth = 6 THz)`；`run_time = 3 ps` |
| 模式规格 | `ModeSpec(num_modes = 2, target_neff = 2.0)`；**Tidy3D 在该 450×220 条上 mode0 = 物理 TE0、mode1 = 物理 TM0**（已按 neff 对照确认） |
| 输出面 | `FieldMonitor`（x-normal）@ x = 33.1 µm，`size = (0, 1.0, 3.0) µm`，`freqs = [c/1.55 µm]`（单频） |
| 输入面 | `FieldMonitor`（x-normal）@ x = -1.5 µm（用于归一化 P_in） |
| 分析量 | 每臂 z 窗口 **±0.45 µm** 内积分 **S_x = ½Re(E_y H_z\* − E_z H_y\*)** + 该窗口 **E_x / E_y 能量占比**（偏振判别） |

## 6. 波长与成本

| 项 | 值 |
|---|---|
| 频点 | **单频 1.55 µm**（`tidy3d_FIELD_TE.csv` 即该频点结果） |
| 成本 | 单次 3D FDTD 场积分跑（35.5 × 2.8 × 6.0 µm、dl 10 nm 量级）**按 `estimate_cost` 报价后再提交**；本报告数据来自已完成的运行 |

## 7. 已知事项（复现时请留意）

1. **必须用场积分法**：2 模式端口法在强辐射/多模耦合结构上会给出伪影（TM 通道曾出现"两臂均 ≈0"的假结果）；
   场积分法不依赖端口模式序号。
2. **2D FDE 与 3D 的 L_c 差异大**：FDE 给 L_c(TE) = 15.8 µm，3D 数据反推 **L_c(TE) ≈ 130 µm**（≫ L = 31.5 µm）
   ⇒ 最终 L / gap 必须在 3D 内重定标后才谈达标。
3. **模式标签反号**：Tidy3D 的 mode0/mode1 与 Lumerical 的标签相反；一律**按 neff** 认模（本项目规范）。
4. **参考直波导对照必须一起跑**：λ = 1503 nm 处参考直波导自身给出 neff = 1.4891（近包层）⇒ 该频点的数据不可信；
   仅 λ ≥ 1519 nm 的对照通过（TE0 99.7%、TM0 92.0%）。
5. 长任务请以后台进程运行并轮询日志；stdout 需 `python -u`。

## 8. 数据依赖（自洽性）

| 项 | 位置 | 说明 |
|---|---|---|
| **离线重分析输入** | **`../data/SPLIT_FIELD_TE.hdf5`（801,440 B ≈ 0.78 MB，已纳入本仓）** | Tidy3D 场监视器原始数据；`--hdf5` 指向它即可 **0 credit** 复算 §5 的全部分析量 |
| 原始运行目录 | 本机项目临时目录（未纳入本仓） | 仓内副本与它 **SHA-256 前 16 位一致**（`70efca91ddf3652d`）⇒ **不依赖仓外文件也能复现** |
| 仓外不可纳入的依赖 | **无** | 本器件的分析输入均已随仓提供（CSV + 场图 + hdf5） |

> 维护提醒：脚本里的 `path=os.path.join(OUT, name + ".hdf5")` 是**新跑一次**时的写出路径（落在你本机的临时目录）；
> **离线重分析**请用上表里仓内的 hdf5。
