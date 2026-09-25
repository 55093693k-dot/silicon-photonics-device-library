# -*- coding: utf-8 -*-
"""run_bilevel_psr.py —— 复刻公开设计 **BilevelPSR**（taper/绝热耦合路线，C 波段）

公开来源：https://docs.flexcompute.com/projects/tidy3d/en/v2.11.0/notebooks/BilevelPSR.html
几何 / 材料 / 网格 / 源 / 端口的完整定义都写在本脚本内（见下方 P 与 make_* 函数），
设置说明见同目录 `SIMULATION_SETUP.md`，结论与引用限制见上一级 `README.md`。

阶段：
  --stage geometry   出结构图（俯视 + 横截面，含尺寸标注）—— 跑仿真前先核对几何
  --stage modes      本地 ModeSolver 在该示例的 4 个 x 位置复算模式，与该示例自带的核对表比对（含 neff>1.44 物理窗口检查）
"""
from __future__ import annotations

import argparse
import os
import sys

import gdstk
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import tidy3d as td
import tidy3d.web as web

for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:                                              # noqa: BLE001
        pass

HERE = os.path.dirname(os.path.abspath(__file__))

# ---------------------------------------------------------------------------
# 该示例参数（逐项与 notebook 一致）
# ---------------------------------------------------------------------------
P = dict(
    L_blt=100.0, L_ac=300.0, L_s=5.0, L_t=30.0,
    w_1=0.45, w_2=0.55, w_pes=1.55, w_3=0.85,
    w_4=0.2, w_5=0.65, w_6=0.5, gap=0.2,
    t_pes=0.09, t_si=0.22,
    R=300.0, theta=np.pi / 30.0,
    sidewall_angle=10 * np.pi / 180.0,
    inf_eff=1e5,
    lda0=1.55, band=(1.50, 1.58), n_pts=9,
)


def si_medium():
    return td.material_library["cSi"]["Palik_Lossless"]


def freqs_of(p):
    return td.C_0 / np.linspace(p["band"][0], p["band"][1], p["n_pts"])


def make_cell(p):
    """该示例的 gdstk 定义：上波导 + 下波导（逐行照抄 notebook）。"""
    inf_eff = p["inf_eff"]
    L_blt, L_ac, L_s, L_t = p["L_blt"], p["L_ac"], p["L_s"], p["L_t"]
    w_1, w_2, w_3, w_5 = p["w_1"], p["w_2"], p["w_3"], p["w_5"]
    w_4, w_6, gap = p["w_4"], p["w_6"], p["gap"]
    R, theta = p["R"], p["theta"]

    cell = gdstk.Cell("device")
    top_wg = gdstk.RobustPath((-inf_eff, 0), w_1, layer=1, datatype=0)
    top_wg.horizontal(0)
    top_wg.horizontal(L_blt / 2, w_2)
    top_wg.horizontal(L_blt, w_3)
    top_wg.horizontal(L_blt + L_s)
    top_wg.segment((L_blt + L_s + L_ac, (w_5 - w_3) / 2), w_5)
    top_wg.horizontal(L_blt + L_s + L_ac + L_t)
    top_wg.arc(R, -np.pi / 2, -np.pi / 2 + theta)
    top_wg.arc(R, np.pi / 2 + theta, np.pi / 2)
    top_wg.horizontal(inf_eff)
    cell.add(top_wg)

    bottom_wg = gdstk.RobustPath((L_blt + L_s, (-w_4 - w_3 - 2 * gap) / 2), w_4,
                                 layer=1, datatype=0)
    bottom_wg.segment((L_blt + L_s + L_ac, (-w_3 - 2 * gap - w_6) / 2), w_6)
    bottom_wg.arc(R, np.pi / 2, np.pi / 2 - theta)
    bottom_wg.arc(R, -np.pi / 2 - theta, -np.pi / 2)
    bottom_wg.horizontal(inf_eff)
    cell.add(bottom_wg)
    return cell


def make_structures(p):
    si = si_medium()
    v = [(0, p["w_1"] / 2), (p["L_blt"] / 2, p["w_pes"] / 2), (p["L_blt"], p["w_3"] / 2),
         (p["L_blt"], -p["w_3"] / 2), (p["L_blt"] / 2, -p["w_pes"] / 2), (0, -p["w_1"] / 2)]
    pes = td.Structure(
        geometry=td.PolySlab(vertices=v, axis=2, slab_bounds=(0, p["t_pes"]),
                             sidewall_angle=p["sidewall_angle"]),
        medium=si, name="partially_etched_slab")
    geos = td.PolySlab.from_gds(make_cell(p), gds_layer=1, axis=2,
                                 slab_bounds=(0, p["t_si"]),
                                 sidewall_angle=p["sidewall_angle"])
    wgs = [td.Structure(geometry=g, medium=si, name=f"wg_{i}") for i, g in enumerate(geos)]
    return [pes] + wgs


# ---------------------------------------------------------------------------
# 结构图
# ---------------------------------------------------------------------------
def draw_xsec_mid(ax, p):
    """x=L_blt/2 截面：部分刻蚀平板（w_pes @ t_pes）+ 上波导（w_2 @ t_si）。"""
    ax.add_patch(plt.Rectangle((-1.6, -0.5), 3.2, 1.1, fc="#cfe3f5", ec="none", zorder=0))
    ax.add_patch(plt.Rectangle((-p["w_pes"] / 2, 0), p["w_pes"], p["t_pes"],
                               fc="#e6a17b", ec="k", lw=0.5, zorder=1))
    ax.add_patch(plt.Rectangle((-p["w_2"] / 2, 0), p["w_2"], p["t_si"],
                               fc="#c94f2b", ec="k", lw=0.6, zorder=2))
    dim_h(ax, -p["w_pes"] / 2, p["w_pes"] / 2, p["t_si"] + 0.12,
          f"w_pes={p['w_pes']*1000:.0f}nm (partial etch slab)")
    dim_h(ax, -p["w_2"] / 2, p["w_2"] / 2, -0.20, f"w_2={p['w_2']*1000:.0f}nm")
    dim_v(ax, 0, p["t_si"], 1.30, f"t_si={p['t_si']*1000:.0f}nm", dx=-0.35)
    dim_v(ax, 0, p["t_pes"], 1.05, f"t_pes={p['t_pes']*1000:.0f}nm", dx=-0.30)
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-0.55, 0.8)
    ax.set_xlabel("y (µm)")
    ax.set_ylabel("z (µm)")
    ax.set_title(f"Cross-section @ x=L_blt/2 = {p['L_blt']/2:.0f} µm")


def draw_xsec_coupler(ax, p):
    """x=L_blt+L_s+L_ac 截面：上波导 w_5（中心 y=(w_5-w_3)/2）+ 下波导 w_6，间隙 gap。"""
    y_top = (p["w_5"] - p["w_3"]) / 2
    y_bot = (-p["w_3"] - 2 * p["gap"] - p["w_6"]) / 2
    ax.add_patch(plt.Rectangle((-2.0, -1.6), 4.0, 2.0, fc="#cfe3f5", ec="none", zorder=0))
    ax.add_patch(plt.Rectangle((y_top - p["w_5"] / 2, 0), p["w_5"], p["t_si"],
                               fc="#c94f2b", ec="k", lw=0.6, zorder=2))
    ax.add_patch(plt.Rectangle((y_bot - p["w_6"] / 2, 0), p["w_6"], p["t_si"],
                               fc="#c94f2b", ec="k", lw=0.6, zorder=2))
    dim_h(ax, y_top - p["w_5"] / 2, y_top + p["w_5"] / 2, p["t_si"] + 0.10,
          f"w_5={p['w_5']*1000:.0f}nm (top)")
    dim_h(ax, y_bot - p["w_6"] / 2, y_bot + p["w_6"] / 2, -0.18, f"w_6={p['w_6']*1000:.0f}nm (bottom)")
    dim_h(ax, y_top - p["w_5"] / 2, y_bot + p["w_6"] / 2, 0.42,
          f"gap={p['gap']*1000:.0f}nm  (edge-to-edge)", dy=0.02)
    dim_v(ax, 0, p["t_si"], -1.85, f"t_si={p['t_si']*1000:.0f}nm", dx=-0.25)
    ax.set_xlim(-2.0, 2.0)
    ax.set_ylim(-1.6, 0.9)
    ax.set_xlabel("y (µm)")
    ax.set_ylabel("z (µm)")
    ax.set_title(f"Cross-section @ end of adiabatic coupler (x={p['L_blt']+p['L_s']+p['L_ac']:.0f} µm)")


def dim_h(ax, x0, x1, y, text, fs=7.5, dy=0.0):
    from matplotlib.patches import FancyArrowPatch
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="<->", mutation_scale=8,
                                 color="#333333", lw=1.0, shrinkA=0, shrinkB=0))
    ax.text((x0 + x1) / 2, y + dy, text, ha="center", va="bottom", fontsize=fs)


def dim_v(ax, y0, y1, x, text, fs=7.5, dx=0.0):
    from matplotlib.patches import FancyArrowPatch
    ax.add_patch(FancyArrowPatch((x, y0), (x, y1), arrowstyle="<->", mutation_scale=8,
                                 color="#333333", lw=1.0, shrinkA=0, shrinkB=0))
    ax.text(x + dx, (y0 + y1) / 2, text, ha="left", va="center", fontsize=fs)


def small_sim(p, x_max=None, y_span=7.0, z_span=1.75):
    """用于绘图/模式求解的**局部**仿真（全长 435 µm 的全域不在此处构建：见 NOTES 的 EME 说明）。"""
    si = si_medium()
    _, sio2 = None, td.material_library["SiO2"]["Palik_Lossless"]
    x_max = x_max if x_max is not None else (p["L_blt"] + p["L_s"] + p["L_ac"] + p["L_t"] + 5)
    return td.Simulation(
        center=((x_max - 5) / 2, 0, 0.11), size=(x_max + 5, y_span, z_span),
        grid_spec=td.GridSpec.auto(min_steps_per_wvl=20, wavelength=p["lda0"]),
        structures=make_structures(p), run_time=1e-12,
        boundary_spec=td.BoundarySpec.all_sides(boundary=td.PML()), medium=sio2)


def full_sim(p):
    """**整器件**仿真域（仅用于结构绘图；不构建网格，故 ~500 µm 也无成本）。
    y 跨度取 **±6 µm**（12 µm），避免边缘波导被截断。"""
    sio2 = td.material_library["SiO2"]["Palik_Lossless"]
    x_end = p["L_blt"] + p["L_s"] + p["L_ac"] + p["L_t"] + 2 * p["R"] * p["theta"] + 20
    return td.Simulation(
        center=(x_end / 2 - 5, 0, 0.11), size=(x_end + 10, 12.0, 1.75),
        grid_spec=td.GridSpec.auto(min_steps_per_wvl=20, wavelength=p["lda0"]),
        structures=make_structures(p), run_time=1e-12,
        boundary_spec=td.BoundarySpec.all_sides(boundary=td.PML()), medium=sio2), x_end


def geometry_bounds(p):
    """打印结构总包围盒（用于确认 y 向是否被截断）。"""
    from tidy3d import GeometryGroup
    geos = [s.geometry for s in make_structures(p)]
    try:
        b = GeometryGroup(geometries=geos).bounds
        print("structure bounds (xmin,xmax,ymin,ymax,zmin,zmax) =",
              tuple(round(v, 3) for v in (b[0][0], b[1][0], b[0][1], b[1][1], b[0][2], b[1][2])))
    except Exception as e:                                          # noqa: BLE001
        print("bounds failed:", repr(e))


def stage_geometry(p):
    """三张图：① 整器件俯视（两层）② 分段放大俯视 ③ 横截面（带尺寸）。"""
    figs = os.path.join(HERE, "figs")
    os.makedirs(figs, exist_ok=True)
    sim, x_end = full_sim(p)
    x_lo, x_hi = -5.0, x_end
    geometry_bounds(p)          # 打印结构包围盒，核对 y 域 ±6 µm 是否完整覆盖

    # ① 整器件俯视（x 与 y 比例不同，图注已说明）
    fig, axs = plt.subplots(2, 1, figsize=(24, 7), dpi=130)
    for ax, z, tag in ((axs[0], p["t_pes"] / 2, "partial-etch slab layer"),
                       (axs[1], p["t_si"] / 2, "silicon waveguide layer")):
        sim.plot(z=z, ax=ax)
        ax.set_aspect("auto")
        ax.set_xlim(x_lo, x_hi)
        ax.set_title(f"WHOLE DEVICE top view — z-slice at z={z:.3f} µm ({tag}). "
                     f"NOTE: single Si layer (t_si=220 nm) + 90-nm partial etch; "
                     f"this is a SLICE, not a second thickness. "
                     f"Total length ≈ {x_hi - x_lo:.0f} µm (x:y aspect not to scale)")
    f1 = os.path.join(figs, "bilevel_psr_topview_FULL.png")
    fig.tight_layout()
    fig.savefig(f1, bbox_inches="tight")
    plt.close(fig)

    # ② 分段放大俯视（同一层，4 段覆盖整器件）
    segs = [(x_lo, 110, "1) input + bi-level taper (w_1=0.45 → w_2=0.55 → w_3=0.85 µm)"),
            (100, 215, "2) straight (L_s=5) + adiabatic coupler entry (bottom wg w_4=0.2 µm)"),
            (205, 410, "3) adiabatic coupler (upper w_5=0.65 / lower w_6=0.5 / gap=0.2 µm)"),
            (400, x_hi, "4) output S-bends (R=300 µm) and output waveguides")]
    fig, axs = plt.subplots(len(segs), 1, figsize=(16, 16), dpi=130)
    for ax, (x0, x1, name) in zip(axs, segs):
        sim.plot(z=p["t_si"] / 2, ax=ax)
        ax.set_xlim(x0, x1)
        ax.set_aspect("auto")
        ax.set_title(name, fontsize=10)
    f2 = os.path.join(figs, "bilevel_psr_topview_ZOOM.png")
    fig.tight_layout()
    fig.savefig(f2, bbox_inches="tight")
    plt.close(fig)

    # ③ 横截面（带尺寸标注）
    fig, axs = plt.subplots(1, 2, figsize=(13, 5), dpi=130)
    draw_xsec_mid(axs[0], p)
    draw_xsec_coupler(axs[1], p)
    f3 = os.path.join(figs, "bilevel_psr_xsec.png")
    fig.tight_layout()
    fig.savefig(f3, bbox_inches="tight")
    plt.close(fig)

    for f in (f1, f2, f3):
        print("figure ->", f, os.path.getsize(f), "bytes")
    return f1, f2, f3


# ---------------------------------------------------------------------------
# 本地模式复核（免费）：在该示例的 4 个 x 位置复算模式，与该示例自带的核对表比对
# ---------------------------------------------------------------------------
def xsec_model(p, x):
    """按该示例参数给出 x 处的**局部截面**（宽度随 x 线性/分段插值）。"""
    y = {}
    if x <= p["L_blt"]:
        f = x / (p["L_blt"] / 2) if x <= p["L_blt"] / 2 else 1 + (x - p["L_blt"] / 2) / (p["L_blt"] / 2)
        y["w_slab"] = (p["w_1"] + (p["w_pes"] - p["w_1"]) * f) if x <= p["L_blt"] / 2 \
            else (p["w_pes"] + (p["w_3"] - p["w_pes"]) * (f - 1))
    if x <= p["L_blt"] / 2:
        y["w_top"] = p["w_1"] + (p["w_2"] - p["w_1"]) * (x / (p["L_blt"] / 2))
    elif x <= p["L_blt"]:
        y["w_top"] = p["w_2"] + (p["w_3"] - p["w_2"]) * ((x - p["L_blt"] / 2) / (p["L_blt"] / 2))
    elif x <= p["L_blt"] + p["L_s"]:
        y["w_top"] = p["w_3"]
    else:
        f = (x - (p["L_blt"] + p["L_s"])) / p["L_ac"]
        y["w_top"] = p["w_3"] + (p["w_5"] - p["w_3"]) * f
        y["y_top"] = (p["w_5"] - p["w_3"]) / 2 * f
    if x >= p["L_blt"] + p["L_s"]:
        f = (x - (p["L_blt"] + p["L_s"])) / p["L_ac"]
        y["w_bot"] = p["w_4"] + (p["w_6"] - p["w_4"]) * f
        y["y_bot"] = ((-p["w_4"] - p["w_3"] - 2 * p["gap"]) / 2
                      + ((-p["w_3"] - 2 * p["gap"] - p["w_6"]) / 2
                         - (-p["w_4"] - p["w_3"] - 2 * p["gap"]) / 2) * f)
    return y


def xsec_sim_plane(p, x, y_c=0.0):
    """**Tidy3D 约定**：器件沿 x 传播 → 模式平面为 **x-normal**（面内轴 = y(横向分离), z(高度)）。
    局部截面按该示例的局部尺寸构建；仿真域 y 取 **±6 µm**，平面取 ±1.5 µm。"""
    si = si_medium()
    sio2 = td.material_library["SiO2"]["Palik_Lossless"]
    m = xsec_model(p, x)
    structs = []
    if "w_slab" in m:
        structs.append(td.Structure(geometry=td.Box(center=(0, y_c, p["t_pes"] / 2),
                                                    size=(3, m["w_slab"], p["t_pes"])), medium=si))
    structs.append(td.Structure(geometry=td.Box(center=(0, y_c + m.get("y_top", 0.0), p["t_si"] / 2),
                                                size=(3, m["w_top"], p["t_si"])), medium=si))
    if "w_bot" in m:
        structs.append(td.Structure(geometry=td.Box(center=(0, y_c + m["y_bot"], p["t_si"] / 2),
                                                    size=(3, m["w_bot"], p["t_si"])), medium=si))
    sim = td.Simulation(center=(0, 0, 0), size=(3.0, 12.0, 1.75),
                        grid_spec=td.GridSpec.auto(min_steps_per_wvl=20, wavelength=p["lda0"]),
                        structures=structs, run_time=1e-13,
                        boundary_spec=td.BoundarySpec.all_sides(boundary=td.PML()), medium=sio2)
    plane = td.Box(center=(0, 0, p["t_si"] / 2), size=(0.0, 3.0, 1.4))   # x-normal（Tidy3D 约定）
    return sim, plane


def stage_modes(p):
    """该示例自带的核对表：x=0 → TE0/TM0；L_blt/2 → TE0/混合；L_blt+L_s → TE0/TE1；耦合器末端 → 上下波导 TE0。"""
    positions = [0.0, p["L_blt"] / 2, p["L_blt"] + p["L_s"],
                 p["L_blt"] + p["L_s"] + p["L_ac"]]
    from tidy3d.plugins.mode import ModeSolver
    for x in positions:
        sim, plane = xsec_sim_plane(p, x)
        spec = td.ModeSpec(num_modes=4, target_neff=2.6).updated_copy(num_pml=(6, 6))
        data = ModeSolver(simulation=sim, plane=plane, mode_spec=spec,
                          freqs=[td.C_0 / p["lda0"]]).solve()
        neff = np.asarray(data.n_eff.values).ravel()
        pf = data.pol_fraction
        te = tm = None
        try:
            te = np.asarray(pf.sel(pol="te").values).ravel()
            tm = np.asarray(pf.sel(pol="tm").values).ravel()
        except Exception as e:                                     # noqa: BLE001
            print(f"  !! pol_fraction 读取失败（易错项）：{e!r}")
            print(f"     dims={getattr(pf, 'dims', None)}  coords={list(getattr(pf, 'coords', {}))}")
            arr = np.asarray(pf.values)
            if arr.ndim >= 2 and arr.shape[-1] == 2:               # 兜底：最后一轴 = (te, tm)
                te, tm = arr[..., 0].ravel(), arr[..., 1].ravel()
                print("     → 已用兜底方式（最后一轴=(te,tm)）读取；请核对该假设")
            else:
                print("     → 兜底失败：**本次模式识别仅按 neff，偏振过滤未生效**（必须在此处停下修正）")
        print(f"--- x = {x:.1f} µm ---")
        for i, n in enumerate(neff):
            flag = "OK" if n > 1.44 else "!! neff<n_clad(1.44) → 泄漏模（易错项）"
            extra = "" if te is None else f"  te={te[i]:.3f} tm={tm[i]:.3f}"
            print(f"  mode{i}: neff={n:.4f}{extra}   {flag}")


def eme_explicit_grid(p, num_modes=4):
    """第二次修正（**唯一单一变量 = EME 网格**；其余设置与上一轮完全一致）。

    依据 Tidy3D API（docs.flexcompute.com, v2.9.3/v2.11.2）：
      `EMEExplicitGrid(mode_specs=[每格一个], boundaries=[内部边界; 个数 = 格数-1; 严格递增])`
    分区（x 从 -10 µm 起）：
      ① bi-level taper 0 → L_blt：2 µm/格（profile 缓变）
      ② 突变界面（L_blt−5 → L_blt+L_s+5）：0.5 µm/格（90 nm 台阶 + 下波导骤现 w_4=0.2 µm）
      ③ 绝热耦合器 L_blt+L_s → x_dev：5 µm/格
    """
    x_dev = p["L_blt"] + p["L_s"] + p["L_ac"] + p["L_t"] + 2 * p["R"] * p["theta"]
    x_start, x_end_sim = -10.0, x_dev + 10.0
    x_tap = p["L_blt"]
    x_junc = p["L_blt"] + p["L_s"]

    xs = []
    xs += list(np.arange(0.0, x_tap + 1e-6, 2.0))                     # ① 2 µm
    xs += list(np.arange(x_tap - 5.0, x_junc + 5.0 + 1e-6, 0.5))      # ② 0.5 µm
    xs += list(np.arange(x_junc, x_dev + 1e-6, 5.0))                  # ③ 5 µm
    bnd = sorted({round(float(v), 5) for v in xs
                  if x_start + 1e-6 < float(v) < x_end_sim - 1e-6})   # 只保留内部边界
    specs = [td.EMEModeSpec(num_modes=num_modes) for _ in range(len(bnd) + 1)]
    print(f"EMEExplicitGrid: {len(bnd) + 1} cells  ({len(bnd)} internal boundaries) "
          f"| taper 2 µm / junction 0.5 µm / coupler 5 µm")
    return td.EMEExplicitGrid(mode_specs=specs, boundaries=tuple(bnd))


def make_eme_sim(p, num_cells=60, num_modes=4):
    """Tidy3D EME 配方：axis=0（x 向传播）+ EME 网格 + port_offsets + 系数监视器。
    此后 eme_grid_spec 改为 `EMEExplicitGrid`（在几何突变处分界），不再是 EMEUniformGrid(60)。"""
    sio2 = td.material_library["SiO2"]["Palik_Lossless"]
    x_dev = p["L_blt"] + p["L_s"] + p["L_ac"] + p["L_t"] + 2 * p["R"] * p["theta"]
    pad = 20.0
    return td.EMESimulation(
        center=(x_dev / 2, 0, 0.11), size=(x_dev + 2 * pad, 12.0, 1.75),   # y=±6 µm
        medium=sio2, structures=make_structures(p),
        axis=0, freqs=freqs_of(p),
        grid_spec=td.GridSpec.auto(min_steps_per_wvl=20, wavelength=p["lda0"]),
        eme_grid_spec=eme_explicit_grid(p, num_modes=num_modes),   # 显式网格（在几何突变处分界）
        port_offsets=(pad / 2, pad / 2),
        constraint="unitary",          # Tidy3D 建议：抑制不稳定辐射模
        store_port_modes=True,         # 关键诊断：保存端口模式（用于判定 mode_index ↔ 波导）
        monitors=[td.EMEModeSolverMonitor(name="port_modes", num_modes=num_modes,
                                         size=(0, td.inf, td.inf), center=(x_dev / 2, 0, 0.11)),
                  td.EMEModeSolverMonitor(name="modes_in", num_modes=num_modes,
                                          size=(0, td.inf, td.inf), center=(-pad / 2, 0, 0.11)),
                  td.EMEModeSolverMonitor(name="modes_out", num_modes=num_modes,
                                          size=(0, td.inf, td.inf), center=(x_dev + pad / 2, 0, 0.11)),
                  td.EMECoefficientMonitor(name="coeffs", size=(td.inf, td.inf, td.inf))],
    )


def stage_eme(p, submit=False):
    """EME 验证：输入 TE0/TM0 注入 → 读输出端两个波导的 TE0/TM0 透射（判 TM→TE>90%）。"""
    sim = make_eme_sim(p, num_modes=2)   # 单一变量：端口模式 4→2（去掉近截止/频率跳变模）
    print("=== 几何核对报告 ===")
    print(" 传播方向 = +x | 分离方向 = y（两波导横向分开）| 高度方向 = z（Si 220nm + 90nm 部分刻蚀，单层）")
    print(f" 端口平面 = **x-normal**（垂直于 +x）| 仿真域 size={sim.size} center={sim.center}")
    print(f" 模式索引：输入 mode_index_in 0=TE0 / 1=TM0（按 neff 识别）；输出端按 neff 区分上下波导")
    print("=== 提交前校验 ===")
    sim.validate_pre_upload()
    print(" validate_pre_upload: PASS")
    task = f"bilevel_psr_eme_Lac{p['L_ac']:.0f}"
    job = web.Job(simulation=sim, task_name=task, verbose=False)
    print(" estimate_cost (FlexCredit):", web.estimate_cost(job.task_id))
    if not submit:
        print(" [dry-run] 未提交；加 --submit 才提交（会先打印成本估算）。")
        return None
    data = job.run(path=os.path.join(HERE, f"data_{task}.hdf5"))
    print(" === smatrix 结构 ===")
    try:
        S = data.smatrix
        print(" smatrix dims:", list(getattr(S, "dims", {})))
        print(" S21 dims:", list(getattr(S.S21, "dims", {})), " shape:", getattr(S.S21, "shape", None))
        print(" coords:", {k: (list(v.values)[:6] if hasattr(v, "values") else v)
                           for k, v in getattr(S.S21, "coords", {}).items()})
    except Exception as e:                                          # noqa: BLE001
        print(" smatrix 读取失败：", repr(e))
    try:
        pm = data["port_modes"]
        print(" port_modes neff (看能否区分上下波导):", np.asarray(pm.n_eff.values))
    except Exception as e:                                          # noqa: BLE001
        print(" port_modes 读取失败：", repr(e))
    print(" data file:", os.path.join(HERE, f"data_{task}.hdf5"))
    return data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="geometry", choices=["geometry", "modes", "eme"])
    ap.add_argument("--submit", action="store_true", help="EME：真正提交（默认 dry-run 只报成本）")
    a = ap.parse_args()
    print("BilevelPSR 复刻参数：", {k: round(v, 4) for k, v in P.items()
                                    if isinstance(v, float)})
    if a.stage == "geometry":
        stage_geometry(P)
    elif a.stage == "modes":
        stage_modes(P)
    else:
        stage_eme(P, submit=a.submit)


if __name__ == "__main__":
    main()


