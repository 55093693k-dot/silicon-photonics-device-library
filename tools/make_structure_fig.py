# -*- coding: utf-8 -*-
"""make_structure_fig.py —— 结构示意图生成器（俯视图 + 横截面，含尺寸标注）

用途：当某个器件只有几何参数、没有现成的结构图时，用这份脚本按参数画出
"俯视图 + 横截面"两张带尺寸标注的图。**所有图必须带尺寸标注（µm / nm）**。

证据等级提醒：这种图是**按几何参数重绘的示意图**，不是求解器的结构视图。
如果该器件已经有**真实仿真场图 / 模式图**，优先直接使用真实图（证据更强），
本脚本只在没有真实图时使用。

用法：
  python make_structure_fig.py --list
  python make_structure_fig.py --device splitter --out ../psr/split_TE_TM_DC
  python make_structure_fig.py --device delay_line --out ../delay_line
"""
from __future__ import annotations

import argparse
import math
import os
import sys

# console code pages are not always UTF-8; never let a print() kill the script
for _s in (sys.stdout, sys.stderr):
    try:
        _s.reconfigure(encoding="utf-8", errors="replace")
    except Exception:                                                # noqa: BLE001
        pass

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Polygon, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
try:
    import fig_channel as FC          # optional: adds the figure-channel label
except Exception:                     # noqa: BLE001  (script still works without it)
    FC = None

SI = "#c94f2b"      # 硅
OX = "#cfe3f5"      # SiO2

# ---- 器件几何（µm，与各器件 README 的参数表一致）-------------------------
DEVICES = {
    # PSR 分束器（见 psr/split_TE_TM_DC/README.md §2）
    "splitter": dict(id="PSR_SPLIT_TE_TM_DC_v1", h=0.220, w1=0.450, w2=0.450,
                     gap=0.150, L=31.5, stub=2.0, ramp=3.0, off=1.35),
    # 1 m 延迟线（见 delay_line/README.md）
    "delay_line": dict(id="DELAY_1M_SOI_v1", h=0.220, w=0.450, n_turn=6,
                       pitch=40.0, span=600.0),
    # 绝热双台阶锥偏振旋转器（官方 BilevelPSR 适配版；
    # 见 psr/rotate_adiabatic_BLT/README.md §2 —— 图中不含任何未认证性能数值）
    "psr_rotator": dict(id="PSR_ROTATE_ADIABATIC_BLT_v1",
                        t_si=0.220, t_pes=0.090, w_1=0.450, w_2=0.550, w_3=0.850,
                        w_4=0.200, w_5=0.650, w_6=0.500, gap=0.200, w_pes=1.550,
                        L_blt=100.0, L_s=5.0, L_ac=300.0, L_t=30.0, R=300.0, th_deg=6.0),
}


def dim_h(ax, x0, x1, y, text, fs=8, dy=0.0, color="#333333"):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="<->", mutation_scale=8,
                                 color=color, lw=1.0, shrinkA=0, shrinkB=0))
    ax.text((x0 + x1) / 2, y + dy, text, ha="center", va="bottom", fontsize=fs, color=color)


def dim_v(ax, y0, y1, x, text, fs=8, dx=0.0, color="#333333"):
    ax.add_patch(FancyArrowPatch((x, y0), (x, y1), arrowstyle="<->", mutation_scale=8,
                                 color=color, lw=1.0, shrinkA=0, shrinkB=0))
    ax.text(x + dx, (y0 + y1) / 2, text, ha="left", va="center", fontsize=fs, color=color)


def save_panel(drawer, path, figsize=(9, 5), provenance=""):
    """Draw one panel and save it.

    These are **matplotlib re-draws of the model geometry** -- i.e. the
    *external* channel, never a solver export. If `fig_channel` is available the
    PNG therefore carries the channel label (on-figure footer + PNG tEXt chunk);
    without it the figure is still produced, just without the machine-readable
    label.
    """
    fig, ax = plt.subplots(figsize=figsize, dpi=140)
    drawer(ax)
    fig.tight_layout()
    if FC is not None:
        FC.save(fig, path, "external",
                provenance or "structure figure re-drawn from the model geometry",
                bbox_inches="tight")
    else:
        fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def splitter_top(ax):
    p = DEVICES["splitter"]
    z1 = -(p["gap"] / 2 + p["w1"] / 2)
    z2 = +(p["gap"] / 2 + p["w2"] / 2)
    x0, x1 = -p["stub"], p["L"] + p["stub"]
    ax.add_patch(Rectangle((x0, z1 - 2.2), x1 - x0, (z2 - z1) + 4.4, fc=OX, ec="none", zorder=0))
    ax.add_patch(Rectangle((x0, z1 - p["w1"] / 2), x1 - x0, p["w1"], fc=SI, ec="k", lw=0.4, zorder=2))
    verts = [(x1, z2 + p["w2"] / 2), (x1, z2 - p["w2"] / 2), (p["ramp"], z2 - p["w2"] / 2),
             (0.0, z2 + p["off"] - p["w2"] / 2), (0.0, z2 + p["off"] + p["w2"] / 2),
             (p["ramp"], z2 + p["w2"] / 2)]
    ax.add_patch(Polygon(verts, closed=True, fc=SI, ec="k", lw=0.4, zorder=2))
    ax.text(0.0, z2 + p["off"] + 0.95, "wg2 S-bend in (offset 1.35 µm)", fontsize=8, ha="center")
    ax.text(p["L"] * 0.5, z1 - 1.15, "wg1 straight (bar arm)", fontsize=8, ha="center")
    ax.text(p["L"] * 0.5, z2 + 1.25, "coupling region (w2 = 450 nm)", fontsize=8, ha="center")
    ax.text(p["L"] * 0.65, 0.0, "gap 150 nm", fontsize=7, ha="center", va="center", zorder=3)
    dim_h(ax, 0, p["L"], z2 + 2.5, "L = 31.5 µm", dy=0.06)
    dim_h(ax, -p["stub"], 0, z1 - 2.8, "stub 2.0 µm", dy=-0.9)
    dim_h(ax, 0, p["ramp"], z1 - 4.1, "ramp 3.0 µm", dy=-0.9)
    dim_h(ax, -p["w1"] / 2, p["w1"] / 2, z1 + 0.5, "450 nm", fs=7, dy=1.5)
    dim_h(ax, p["gap"] + p["w2"] / 2, p["gap"] + p["w2"] * 1.5, z2 - 0.5, "450 nm", fs=7, dy=-1.5)
    ax.set_xlim(x0 - 1, x1 + 2)
    ax.set_ylim(z1 - 5.6, z2 + p["off"] + 2.3)
    ax.set_title("Top view (xz) — TE/TM splitter: equal-width DC + length selection")
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("z (µm)")


def splitter_xsec(ax):
    p = DEVICES["splitter"]
    ax.add_patch(Rectangle((-2.2, -0.9), 4.4, 1.55, fc=OX, ec="none", zorder=0))
    ax.add_patch(Rectangle((-p["w1"] / 2, 0.0), p["w1"], p["h"], fc=SI, ec="k", lw=0.5, zorder=2))
    ax.add_patch(Rectangle((p["gap"] + p["w2"] / 2, 0.0), p["w2"], p["h"], fc=SI, ec="k", lw=0.5, zorder=2))
    ax.text(-p["w1"] / 2 - 0.12, p["h"] / 2, "wg1", fontsize=8, ha="right", va="center")
    ax.text(p["gap"] + p["w2"] * 1.5 + 0.12, p["h"] / 2, "wg2", fontsize=8, ha="left", va="center")
    dim_h(ax, -p["w1"] / 2, p["w1"] / 2, p["h"] + 0.45, "w1 = 450 nm")
    dim_h(ax, p["gap"] + p["w2"] / 2, p["gap"] + p["w2"] * 1.5, p["h"] + 0.45, "w2 = 450 nm")
    dim_h(ax, p["w1"] / 2, p["gap"] + p["w2"] / 2, -0.45, "gap = 150 nm", dy=-0.32)
    dim_v(ax, 0, p["h"], -1.5, "H = 220 nm", dx=-0.1)
    ax.text(-2.05, -0.78, "SiO2 cladding", fontsize=7, color="#4a6fa5")
    ax.set_xlim(-2.2, 2.2)
    ax.set_ylim(-1.0, 1.5)
    ax.set_title("Cross-section (z-cut at coupling region) — strip waveguides")
    ax.set_xlabel("z (µm)")
    ax.set_ylabel("y (µm)")


def delay_top(ax):
    p = DEVICES["delay_line"]
    ax.add_patch(Rectangle((-0.05 * p["span"], -26), 1.1 * p["span"],
                           p["pitch"] * p["n_turn"] * 0.55 + 30, fc=OX, ec="none", zorder=0))
    y, xs, w = 0.0, 0.0, 3.0
    step = p["pitch"] * 0.45
    for i in range(p["n_turn"]):
        x_end = p["span"] * (0.30 + 0.55 * (i + 1) / p["n_turn"])
        ax.add_patch(Rectangle((xs, y - w / 2), x_end - xs, w, fc=SI, ec="k", lw=0.4, zorder=2))
        ax.add_patch(Rectangle((x_end - w, y - w / 2), w, step + w, fc=SI, ec="k", lw=0.4, zorder=2))
        y += step
        xs = x_end - w
    ax.text(p["span"] * 0.5, -21, "serpentine routing (schematic, not to scale; drawn 3 µm wide)",
            fontsize=8, ha="center")
    ax.text(p["span"] * 0.5, y + 9,
            "total optical length = 1 m   |   SOI 450 nm strip   |   loss 26.5-28.5 dB",
            fontsize=8, ha="center")
    dim_h(ax, 0, p["span"] * 0.85, -25, "footprint ≈ 600 µm (schematic)", dy=-1.2)
    ax.set_xlim(-0.06 * p["span"], 1.12 * p["span"])
    ax.set_ylim(-31, y + 17)
    ax.set_title("Top view (schematic) — 1 m delay line / routing")
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("z (µm)")


def delay_xsec(ax):
    p = DEVICES["delay_line"]
    ax.add_patch(Rectangle((-1.6, -0.9), 3.2, 1.55, fc=OX, ec="none", zorder=0))
    ax.add_patch(Rectangle((-p["w"] / 2, 0.0), p["w"], p["h"], fc=SI, ec="k", lw=0.5, zorder=2))
    dim_h(ax, -p["w"] / 2, p["w"] / 2, p["h"] + 0.45, "w = 450 nm")
    dim_v(ax, 0, p["h"], -1.05, "H = 220 nm", dx=-0.1)
    ax.text(-1.5, -0.78, "SiO2 cladding", fontsize=7, color="#4a6fa5")
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.0, 1.5)
    ax.set_title("Cross-section — SOI 220 nm strip waveguide")
    ax.set_xlabel("z (µm)")
    ax.set_ylabel("y (µm)")


def _rot_geom(p):
    """器件 x 分段（µm）。几何约定与横截面一致：**两臂间隙恒定 = gap**，
    臂中心 y = ±(gap + w)/2（故臂中心随宽度移动）；x > L_ac 后按 2θ 张开（S-bend 示意）。"""
    x_blt = p["L_blt"]
    x_str = p["L_blt"] + p["L_s"]
    x_ac = x_str + p["L_ac"]
    x_end = x_ac + p["L_t"] + 2 * p["R"] * math.radians(p["th_deg"])
    return x_blt, x_str, x_ac, x_end


def _rot_w_up(x, p):
    """上臂宽度（µm）：0.45 →(0.45 L_blt)→ 0.55 →(0.75 L_blt)→ 0.85 →(L_blt+L_s 后)→ 0.65"""
    x_blt, x_str, _, _ = _rot_geom(p)
    if x <= 0:
        return p["w_1"]
    if x < 0.45 * x_blt:
        return p["w_1"] + (p["w_2"] - p["w_1"]) * x / (0.45 * x_blt)
    if x < 0.75 * x_blt:
        return p["w_2"] + (p["w_3"] - p["w_2"]) * (x - 0.45 * x_blt) / (0.30 * x_blt)
    if x < x_str:
        return p["w_3"]
    return p["w_5"]


def _rot_w_lo(x, p):
    """下臂宽度（µm）：在 L_blt+L_s 处骤现 w_4=0.2，沿 L_ac 线性加宽到 w_6=0.5。"""
    _, x_str, _, _ = _rot_geom(p)
    if x < x_str:
        return 0.0
    f = min(1.0, (x - x_str) / p["L_ac"])
    return p["w_4"] + (p["w_6"] - p["w_4"]) * f


def _rot_pts(p, x0, x1, step):
    """按 step µm 采样 x（限制在器件 x∈[0, x_end] 内）。"""
    _, _, _, x_end = _rot_geom(p)
    xs, x = [], max(0.0, x0)
    while x < min(x_end, x1) - 1e-9:
        xs.append(round(x, 3))
        x += step
    xs.append(round(min(x_end, x1), 3))
    return xs


def _rot_ymed(x, w, sign, p):
    """臂中心 y：±(gap + w)/2，x > L_ac 后按 2θ 张开。"""
    _, _, x_ac, _ = _rot_geom(p)
    y = sign * (p["gap"] + w) / 2
    if x > x_ac:
        y += sign * (x - x_ac) * math.tan(math.radians(p["th_deg"]))
    return y


def rotator_top(ax, seg=None):
    """俯视图。seg=None → **整器件**（x 非等比，长度用尺寸线标出）；
    seg=(x0,x1,cap) → 该段**等比放大**（宽度/间隙可见）。两臂间隙恒为 0.2 µm，保证"分离不重叠"。"""
    p = DEVICES["psr_rotator"]
    x_blt, x_str, x_ac, x_end = _rot_geom(p)
    whole = seg is None
    x0, x1 = (-15.0, x_end + 22.0) if whole else (seg[0], seg[1])
    xs = _rot_pts(p, x0, x1, 2.0 if whole else 0.5)
    # 上臂
    up = [(x, _rot_ymed(x, _rot_w_up(x, p), +1, p) + _rot_w_up(x, p) / 2) for x in xs]
    dn = [(x, _rot_ymed(x, _rot_w_up(x, p), +1, p) - _rot_w_up(x, p) / 2) for x in reversed(xs)]
    # 下臂（x ≥ L_blt+L_s）
    xs2 = [x for x in xs if x >= x_str]
    up2, dn2 = [], []
    if xs2:
        up2 = [(x, _rot_ymed(x, _rot_w_lo(x, p), -1, p) + _rot_w_lo(x, p) / 2) for x in xs2]
        dn2 = [(x, _rot_ymed(x, _rot_w_lo(x, p), -1, p) - _rot_w_lo(x, p) / 2) for x in reversed(xs2)]
    ys = [q[1] for q in up + dn + up2 + dn2]
    ylo_, yhi_ = min(ys), max(ys)
    # 背景
    ax.add_patch(Rectangle((x0, ylo_ - 1.8), x1 - x0, (yhi_ - ylo_) + 3.6, fc=OX, ec="none", zorder=0))
    ax.text(x0 + 0.02 * (x1 - x0), ylo_ - 1.0, "SiO2 cladding", fontsize=7, color="#4a6fa5")
    # 90 nm 部分刻蚀平板（taper→耦合器末端）
    ax.add_patch(Rectangle((0, -p["w_pes"] / 2), min(x_ac, x1) - max(0, x0) if x0 > 0 else x_ac,
                           p["w_pes"], fc="#f3c9a8", ec="none", zorder=1))
    ax.add_patch(Polygon(up + dn, closed=True, fc=SI, ec="k", lw=0.5, zorder=3))
    if up2:
        ax.add_patch(Polygon(up2 + dn2, closed=True, fc=SI, ec="k", lw=0.5, zorder=3))
    # 分段长度尺寸线
    for xa, xb, t in ((0, x_blt, "L_blt = 100 µm"), (x_blt, x_str, "L_s = 5 µm"),
                      (x_str, x_ac, "L_ac = 300 µm"), (x_ac, x_end, "L_t 30 + S-bend")):
        if xb > x0 and xa < x1:
            dim_h(ax, max(xa, x0), min(xb, x1), ylo_ - 0.55, t, dy=0.06, fs=7)
    # 间隙尺寸线（该段含 x_ac 才画）
    if x0 - 1 <= x_ac <= x1 + 1:
        xg = x_ac - 3 if not whole else x_ac + 10
        dim_v(ax, -p["gap"] / 2, p["gap"] / 2, xg, "gap = 200 nm", dx=0.03 * (x1 - x0), fs=7)
    # 臂宽尺寸线（等比视图才看得见）
    if not whole:
        xm = x0 + 0.35 * (x1 - x0)
        wu = _rot_w_up(xm, p)
        dim_h(ax, _rot_ymed(xm, wu, +1, p) - wu / 2, _rot_ymed(xm, wu, +1, p) + wu / 2,
              yhi_ + 0.25, "upper w = %.2f µm" % wu)
        wl = _rot_w_lo(xm, p)
        if wl > 0:
            dim_h(ax, _rot_ymed(xm, wl, -1, p) - wl / 2, _rot_ymed(xm, wl, -1, p) + wl / 2,
                  ylo_ - 0.30, "lower w = %.2f µm" % wl, dy=-0.30)
        Lx, Ly = (x1 - x0), (yhi_ - ylo_)
        mag = max(1.0, round(0.25 * Lx / max(Ly, 1e-6), 1))   # 让图面长宽比 ≈4:1（合规、易读）
        ax.set_title("Top view — x = %.0f … %.0f µm  |  %s" % (x0, x1, seg[2] if len(seg) > 2 else ""),
                     fontsize=9)
        ax.text(0.99, 0.02, "y (width) exaggerated ×%.0f for visibility — all dimension labels are TRUE values"
                % mag, transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5, color="#555555")
        ax.set_ylim(ylo_ - 1.05, yhi_ + 0.85)
        ax.set_aspect(mag)
    else:
        ax.text(0.5 * x_blt, yhi_ + 0.55, "bi-level taper: w_1 0.45 → w_2 0.55 → w_3 0.85 µm",
                fontsize=7, ha="center")
        ax.text(x_str + p["L_ac"] * 0.5, yhi_ + 0.55, "adiabatic coupler (w_5 0.65 / w_6 0.5 µm)",
                fontsize=7, ha="center")
        ax.text(x_str + 2, ylo_ - 0.95, "lower arm appears at x = L_blt + L_s (w_4 = 0.2 µm)", fontsize=7)
        ax.text(x_end + 3, 0.0, "S-bend\nR 300 µm\n2θ 12°", fontsize=7, va="center")
        xcomp = max(1, round((x1 - x0) / max(yhi_ - ylo_, 1e-6) / 4.2))
        ax.text(0.99, 0.02, "x compressed ≈×%d (NOT to scale) — y (widths) as drawn; see z1–z4 for detail views"
                % xcomp, transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5, color="#555555")
        ax.set_title("WHOLE-DEVICE top view — bi-level taper PSR (BilevelPSR replica) | single Si layer "
                     "220 nm + 90 nm partial etch\ntotal length ≈ 537.8 µm | x compressed for readability",
                     fontsize=9)
        ax.set_ylim(ylo_ - 1.35, yhi_ + 1.0)
        ax.set_aspect("auto")
    ax.set_xlim(x0, x1)
    ax.set_xlabel("x (µm)" + ("  [not to scale]" if whole else ""))
    ax.set_ylabel("y (µm)")


def _rot_seg(x0, x1, cap):
    def _draw(ax):
        rotator_top(ax, seg=(x0, x1, cap))
    return _draw


def rotator_xsec(ax):
    """耦合区横截面（后端面）：Si 220 nm + 90 nm 部分刻蚀平板 + 上/下条波导。"""
    p = DEVICES["psr_rotator"]
    ax.add_patch(Rectangle((-1.05, -0.75), 2.1, 1.6, fc=OX, ec="none", zorder=0))
    y_up = (p["gap"] + p["w_5"]) / 2
    y_lo = -(p["gap"] + p["w_6"]) / 2
    ax.add_patch(Rectangle((-p["w_pes"] / 2, 0.0), p["w_pes"], p["t_pes"], fc="#f3c9a8",
                           ec="k", lw=0.4, zorder=1))
    ax.add_patch(Rectangle((y_up - p["w_5"] / 2, 0.0), p["w_5"], p["t_si"], fc=SI, ec="k", lw=0.5, zorder=2))
    ax.add_patch(Rectangle((y_lo - p["w_6"] / 2, 0.0), p["w_6"], p["t_si"], fc=SI, ec="k", lw=0.5, zorder=2))
    dim_h(ax, y_up - p["w_5"] / 2, y_up + p["w_5"] / 2, p["t_si"] + 0.30, "w_5 = 650 nm")
    dim_h(ax, y_lo - p["w_6"] / 2, y_lo + p["w_6"] / 2, p["t_si"] + 0.30, "w_6 = 500 nm")
    dim_h(ax, y_lo + p["w_6"] / 2, y_up - p["w_5"] / 2, p["t_si"] + 0.62, "gap = 200 nm", dy=0.06)
    dim_h(ax, -p["w_pes"] / 2, p["w_pes"] / 2, -0.34, "slab w_pes = 1.55 µm", dy=-0.30, fs=7)
    dim_v(ax, 0, p["t_pes"], p["w_pes"] / 2 + 0.12, "90 nm", dx=0.06, fs=7)
    dim_v(ax, 0, p["t_si"], -p["w_pes"] / 2 - 0.30, "220 nm", dx=0.06, fs=7)
    ax.text(-1.0, -0.70, "SiO2 (n=1.444)", fontsize=7, color="#4a6fa5")
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-0.8, 1.15)
    ax.set_title("Cross-section at the adiabatic coupler (two 220 nm ribs on a 90 nm slab)")
    ax.set_xlabel("lateral y (µm)")
    ax.set_ylabel("height z (µm)")


PANELS = {
    "splitter": [("topview", splitter_top, (10, 6)), ("xsec", splitter_xsec, (8, 4.5))],
    "delay_line": [("topview", delay_top, (10, 6)), ("xsec", delay_xsec, (8, 4.5))],
    # 出图口径（2026-09-11）：只出**整器件**图、固定长方形比例；**不做局部放大**（z1–z4 已移除）
    "psr_rotator": [("topview", rotator_top, (16, 6)), ("xsec", rotator_xsec, (9, 6))],
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--device", choices=sorted(PANELS))
    ap.add_argument("--out", default=".")
    ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    if a.list:
        for k in sorted(PANELS):
            print(k, "->", DEVICES[k]["id"] + "_{topview,xsec}.png")
        return
    if not a.device:
        ap.error("--device required (or --list)")
    os.makedirs(a.out, exist_ok=True)
    dev_id = DEVICES[a.device]["id"]
    for tag, drawer, size in PANELS[a.device]:
        path = os.path.join(a.out, f"{dev_id}_{tag}.png")
        save_panel(drawer, path, size,
                   provenance="structure figure re-drawn from the DEVICES['%s'] "
                              "geometry (not a solver export)" % a.device)
        print("wrote:", path, os.path.getsize(path), "bytes")


if __name__ == "__main__":
    main()
