# 数据清单 — `PSR_SPLIT_TE_TM_DC_v1`

> 用途：让使用者把自己跑出的结果与**本器件的真实数据**逐点对照。
> 数据来源：本项目 MODE FDE（2D 设计筛选）+ Tidy3D 3D FDTD（验证）；数值与 `../README.md` §3 及 `S_PARAMS.csv` 一致。

| 文件 | 内容 | 来源 |
|---|---|---|
| **`S_PARAMS.csv`** | **总表**：① FDE 预测（5 波长）② 3D 参考直波导对照（7 波长，含 neff）③ 3D 场积分结果（1550 nm） | 汇总（各分项见下） |
| `FDE_vs_wavelength.csv` | FDE 的 TE 直通 / TM 交叉 vs 波长原始曲线数据 | Lumerical MODE FDE 超模扫描 |
| `FDE_vs_length.csv` | FDE 的 TE/TM 耦合 vs 耦合长度原始数据（选 L 用） | 同上 |
| `REF_GUIDE_TE.csv` / `REF_GUIDE_TM.csv` | 3D 参考直波导对照：TE0 / TM0 注入的模式幅度与 neff vs 波长 | Tidy3D FDTD（同链路对照） |
| `FIELD_TE.csv` | 3D 场积分结果：各臂能流 S_x、E_x/E_y 占比、功率份额 | Tidy3D FDTD（网格 20/10/10 nm） |
| **`SPLIT_FIELD_TE.hdf5`** | **3D 场监视器原始数据（801,440 B ≈ 0.78 MB，已随仓提供）**——离线重分析输入（**不需访问仓外文件**） | Tidy3D `web.run` 下载产物；SHA-256 前 16 位 `70efca91ddf3652d` |

**场分布图**：`../PSR_SPLIT_TE_TM_DC_v1_mode_maps.png`（TE/TM **超模**真实模式场，FDE 输出）。

**可重跑脚本**：`../sim/`（`run_tidy3d_splitter_field.py` 为主用；设置见 `../sim/SIMULATION_SETUP.md`）。

**口径提醒**：
- `REF_GUIDE_*` 的幅度已换算为**功率**（幅度平方）；λ = 1503 nm 一行是参考直波导**自身失效**的证据，不可用。
- 3D 只有 **λ = 1.55 µm 单频**；FDE 才是 5 波长曲线（2D，仅作设计依据）。
- **离线重分析（0 credit）**：输入 `SPLIT_FIELD_TE.hdf5` 已随仓提供 ⇒ `python ../sim/run_tidy3d_splitter_field.py --hdf5 ../data/SPLIT_FIELD_TE.hdf5 --tag TE` **不依赖仓外文件**即可复算。
- 本器件 **未过 3D 判定线（TE 直通 92.65% < 95%）**，且 **TM 通道与绝对效率未测** ⇒ 任何引用必须附 `../README.md` §5 的适用性限制。
