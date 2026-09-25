# -*- coding: utf-8 -*-
"""analyze_coupler_eme.py -- 离线分析"耦合器段"EME 结果（0 FlexCredit）

目的（补上 run_eme_coupler.py 缺的那一步）：
  1) 从 hdf5 读**端口模式解**（EMEModeSolverMonitor @ x=105 / x=405）→ neff + 场极性
     → **按 neff / 极性认模**（不按位置猜模式身份）；
  2) 打印带标签的 T = |S21|^2（out x in）+ 逐入列和 / 逐出行和；
  3) 输出日志同时落盘（utf-8），避免控制台解码干扰（易错项）。

先用 `python analyze_coupler_eme.py` 跑（不需要 --submit，完全离线）。
"""
import os
import sys
import traceback

import numpy as np
import tidy3d as td

HERE = os.path.dirname(os.path.abspath(__file__))
H5 = os.path.join(HERE, "data_psr_coupler_eme_m4_y6.hdf5")
LOG = os.path.join(HERE, "log_an_cpl.txt")

_LINES = []


def P(s=""):
    print(s)
    _LINES.append(str(s))


def flush():
    with open(LOG, "w", encoding="utf-8") as fh:
        fh.write("\n".join(_LINES) + "\n")


def ycenter(d, comp):
    """每个模式的 |comp|^2 在 y 方向的**重心**（用于判定模式属于上/下波导）。"""
    c = getattr(d, comp, None)
    if c is None or "y" not in c.dims:
        return None
    w = np.abs(np.asarray(c.values)) ** 2
    keep = [i for i, dn in enumerate(c.dims) if dn in ("mode_index", "y")]
    prof = w.sum(axis=tuple(i for i in range(w.ndim) if i not in keep))
    if [c.dims[i] for i in keep][0] == "y":          # 不假设轴序，按维度名判定（易错项）
        prof = np.swapaxes(prof, -1, -2)             # -> (mode, y)
    y = np.asarray(c.coords["y"].values)
    tot = prof.sum(axis=1)
    tot[tot == 0] = 1.0
    return (prof * y[None, :]).sum(axis=1) / tot


def port_table(md, name):
    """打印单个端口的模式表（neff + Ey/Ez/Ex 能量占比），返回 neff 数组。"""
    if name not in md:
        P("[%s] MISSING in monitor_data" % name)
        return None
    d = md[name]
    have = [a for a in ("n_complex", "Ex", "Ey", "Ez") if getattr(d, a, None) is not None]
    P("[%s] available: %s" % (name, have))
    nc = d.n_complex
    P("[%s] n_complex dims=%s shape=%s" % (name, nc.dims, nc.shape))
    neff = np.real(np.asarray(nc.values)).reshape(-1)
    P("[%s] neff = %s" % (name, np.round(neff, 4)))

    frac = {}
    for comp in ("Ey", "Ez", "Ex"):
        c = getattr(d, comp, None)
        if c is None:
            continue
        axes = tuple(i for i, dn in enumerate(c.dims) if dn in ("y", "z"))
        if not axes:
            continue
        frac[comp] = np.abs(np.asarray(c.values)) ** 2        # (f, mode, y, z) 形态
        frac[comp] = frac[comp].sum(axis=axes).reshape(-1)   # -> (f*mode,)
    if frac:
        tot = np.zeros_like(frac[list(frac)[0]], dtype=float)
        for k in frac:
            tot = tot + frac[k]
        tot[tot == 0] = 1.0
        P("[%s] field energy fractions + transverse center of mass (per mode):" % name)
        for i in range(len(tot)):
            row = "  ".join("%s=%.3f" % (k, frac[k][i] / tot[i]) for k in frac)
            cy = ycenter(d, "Ey")
            cz = ycenter(d, "Ez")
            P("   mode %d: %s   neff=%.4f  y_c(Ey)=%s  y_c(Ez)=%s"
              % (i, row, neff[i],
                 ("%.4f" % cy[i]) if cy is not None else "n/a",
                 ("%.4f" % cz[i]) if cz is not None else "n/a"))
    return neff


try:
    P("hdf5: %s" % H5)
    P("exists: %s" % os.path.exists(H5))
    data = None
    for meth in ("from_hdf5", "from_file"):
        f = getattr(td.EMESimulationData, meth, None)
        if f is None:
            continue
        try:
            data = f(H5)
            P("loaded via EMESimulationData.%s()" % meth)
            break
        except Exception as e:                                  # noqa: BLE001
            P("%s failed: %r" % (meth, str(e)[:200]))
    if data is None:
        P("candidates: %s" % [a for a in dir(td.EMESimulationData) if "load" in a or "from" in a])
        flush()
        sys.exit(1)

    md = data.monitor_data
    P("monitor_data keys: %s" % list(md.keys()))

    P("=" * 70)
    neff_in = port_table(md, "modes_in")       # x=105 入口
    P("=" * 70)
    neff_out = port_table(md, "modes_out")     # x=405 出口

    P("=" * 70)
    S = data.smatrix
    S21 = S.S21
    P("S21 dims=%s shape=%s" % (S21.dims, S21.shape))
    order = ("f", "mode_index_out", "mode_index_in")
    # 注意：单频点（f 长度 1）时 squeeze() 会把 f 也删掉 → 用 isel 只去掉 sweep 轴（易错项）
    S21c = S21.isel(sweep_index=0).transpose(*order)
    P("after isel/transpose: dims=%s shape=%s" % (S21c.dims, S21c.shape))
    T = np.abs(np.asarray(S21c.values)) ** 2                   # (nf, out, in)
    lam = td.C_0 / np.asarray(S21c.coords["f"].values).ravel()  # µm
    P("lambda (um): %s" % np.round(lam, 4))
    P("T = |S21|^2  (rows = mode_index_out, cols = mode_index_in):")
    with np.printoptions(precision=5, suppress=False):
        P(str(np.round(T, 5)))
    P("column sums (per INPUT mode, summed over 4 out modes): %s" % np.round(T.sum(axis=1), 5))
    P("row sums    (per OUTPUT mode, summed over 4 in modes): %s" % np.round(T.sum(axis=2), 5))
    P("grand total = %.6f   mean per input column = %.6f"
      % (T.sum(), T.sum() / T.shape[2]))
    if T.shape[0] == 1:
        t = T[0]
        P("KEY conversions (single freq):")
        P("  in0 -> out0 (top TE0  -> top TE0)      = %.5f" % t[0, 0])
        P("  in1 -> out1 (TE1      -> bottom TE0)   = %.5f   <== device function" % t[1, 1])
        P("  in1 -> out0 (TE1      -> top TE0)      = %.5f   (unwanted)" % t[0, 1])
        P("  in0 -> out1 (top TE0  -> bottom TE0)   = %.5f   (crosstalk)" % t[1, 0])

    # 认模（只在拿到端口 neff 时做）
    if neff_in is not None and neff_out is not None:
        P("=" * 70)
        tin = np.real(neff_in).reshape(-1)
        tout = np.real(neff_out).reshape(-1)
        P("input port modes  : %s" % ", ".join("m%d=%.4f" % (i, v) for i, v in enumerate(tin)))
        P("output port modes : %s" % ", ".join("m%d=%.4f" % (i, v) for i, v in enumerate(tout)))
        P("expect @x=105 (from local FDE/FDTD): TE0~2.71, TE1~2.25")
        P("expect @x=405 (from the earlier full-device EME run @1.58um): top-TE0 2.6486, bottom-TE0 2.4976")
    flush()
    P("[done] log written to %s" % LOG)
    flush()
except Exception:                                              # noqa: BLE001
    P("TRACEBACK:")
    P(traceback.format_exc())
    flush()
    raise
