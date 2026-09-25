# tools —— 结构图与性能曲线生成器

两个本地脚本，**0 成本、离线、不联网**。用来在没有现成图件时按几何参数补图。

| 脚本 | 作用 |
|---|---|
| `make_structure_fig.py` | 按几何参数生成"俯视图 + 横截面"两张**带尺寸标注**的结构图（支持 3 个器件） |
| `plot_psr_rotator_perf.py` | 生成偏振旋转器的性能图表（转换效率 vs 波长 + TE0 功率去向） |

依赖：`matplotlib`（`make_structure_fig.py` 只用 stdlib + matplotlib）。

---

## 1. `make_structure_fig.py`

```powershell
python make_structure_fig.py --list
python make_structure_fig.py --device splitter   --out ..\psr\split_TE_TM_DC
python make_structure_fig.py --device delay_line --out ..\delay_line
python make_structure_fig.py --device psr_rotator --out ..\psr\rotate_adiabatic_BLT
```

内置 3 个器件的几何（与各器件 README 的参数表一致）：`splitter`、
`delay_line`、`psr_rotator`。新增器件时在脚本顶部的 `DEVICES` 里加一条即可。

**证据等级提醒（重要）**：这张图是**按几何参数重绘的示意图**，不是求解器的结构视图。
如果该器件已经有**真实仿真场图 / 模式图**，请优先使用真实图 —— 本仓里
`waveguide/`、`coupler/`、`mmi/`、`bend/`、`psr/split_TE_TM_DC/`、`psr/rotate_adiabatic_BLT/`
的图都是**真实图**。本脚本只在没有真实图时用来补位。

**为什么坚持"必须带尺寸标注"**：几何图的作用是让人**读到数**（间隙多少、长度多少），
不带标注的示意图只能看个形状，对复现没有帮助。

---

## 2. `plot_psr_rotator_perf.py`

```powershell
python plot_psr_rotator_perf.py     # 输出 PSR_ROTATE_ADIABATIC_BLT_v1_performance.png 到脚本同目录
```

数据全部写在脚本里（9 个波长的转换效率、TE0 输入功率去向的三个份额）。

⚠️ **这张图的口径限制**（与 `psr/rotate_adiabatic_BLT/README.md` §4 一致）：
它画的是**全器件 EME（`unitary`、2 模式）**的那一组数 —— 那组数在 `passive` 约束下
被证明不可引用，因此这张图**只能看形状趋势**（转换效率随波长的平坦度），
**不能作为绝对性能指标引用**。可引用的数值只有该 README §3.1 的两段：
taper 段 TM0→TE1 = 98.35% 与耦合器段 TE1→下臂 TE0 = 99.233%。
