# -*- coding: utf-8 -*-
"""analyze_coupler_eme_modes.py -- 离线分析任意模式数的"耦合器段"EME 结果（0 FlexCredit）

用法：python analyze_coupler_eme_modes.py --hdf5 <data.hdf5> [--log <out.txt>]

认模规则：**按 neff + 横向场重心**自动定位关键通道，**绝不按模式序号**。
（同目录 `analyze_coupler_eme.py` 是把它写死为 4 模式/固定文件名的版本；本脚本参数化，
供 4 / 6 / 8 / … 模式复用 —— §4 那张 6 模式表就是用本脚本制表的。）

**复用前提（两点，换器件必看）**
1. 本脚本**完全离线**：只用 `tidy3d` 的 `EMESimulationData` 读 hdf5，**不调用 `web`** ⇒
   不需要 API key。运行前提是自备 `tidy3d` 与那一次运行的 `.hdf5`。
2. 下面几处是**本器件标定值**，换器件必须改：
   - `TE1 neff ≈ 2.248`（入口 x = 105 µm 的上波导 TE1）；
   - 出口两臂的**横向场重心** `−0.875`（下臂）/ `−0.10`（上臂）。
   这两个重心来自本器件的几何（w_6 = 0.50 µm / w_5 = 0.65 µm 对应的 y 位置），
   换几何必须重新计算，否则关键通道会被认错。
"""
import argparse
import os
import traceback

import numpy as np
import tidy3d as td

AP = argparse.ArgumentParser()
AP.add_argument("--hdf5", required=True)
AP.add_argument("--log", default=None)
ARGS = AP.parse_args()
H5 = os.path.abspath(ARGS.hdf5)
LOG = os.path.abspath(ARGS.log) if ARGS.log else os.path.splitext(H5)[0] + "_analysis.txt"

_LINES = []


def P(s=""):
    print(s, flush=True)
    _LINES.append(str(s))


def flush():
    with open(LOG, "w", encoding="utf-8") as fh:
        fh.write("\n".join(_LINES) + "\n")


def ycenter(d, comp):
    """每个模式的 |comp|^2 在 y 方向的**重心**（判定上/下波导）。"""
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
    """打印单个端口的模式表（neff + 场能量占比 + 横向重心），返回 neff 数组。"""
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
        frac[comp] = np.abs(np.asarray(c.values)) ** 2        # (f, mode, y, z)
        frac[comp] = frac[comp].sum(axis=axes).reshape(-1)   # -> (f*mode,)
    if frac:
        tot = np.zeros_like(frac[list(frac)[0]], dtype=float)
        for k in frac:
            tot = tot + frac[k]
        tot[tot == 0] = 1.0
        cy = ycenter(d, "Ey")
        cz = ycenter(d, "Ez")
        P("[%s] field energy fractions + transverse center of mass (per mode):" % name)
        for i in range(len(tot)):
            row = "  ".join("%s=%.3f" % (k, frac[k][i] / tot[i]) for k in frac)
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
        raise SystemExit(1)

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
    s = S21
    if "sweep_index" in s.dims:
        s = s.isel(sweep_index=0)
        P("used isel(sweep_index=0)")
    else:
        P("no sweep_index dim -> used as-is (单频点；易错项：不用 squeeze)")
    order = [d for d in ("f", "mode_index_out", "mode_index_in") if d in s.dims]
    s = s.transpose(*order)
    P("after transpose: dims=%s shape=%s" % (s.dims, s.shape))
    T = np.abs(np.asarray(s.values)) ** 2                      # (nf, out, in)
    n_out, n_in = T.shape[1], T.shape[2]
    P("T shape (f,out,in) = %s ; modes/port: in=%d out=%d" % (T.shape, n_in, n_out))
    P("T = |S21|^2  (rows = mode_index_out, cols = mode_index_in):")
    with np.printoptions(precision=5, suppress=False):
        P(str(np.round(T, 5)))
    P("column sums (per INPUT mode, over %d out modes): %s" % (n_out, np.round(T.sum(axis=1), 5)))
    P("row sums    (per OUTPUT mode, over %d in modes): %s" % (n_in, np.round(T.sum(axis=2), 5)))
    P("grand total = %.6f | mean per input column = %.6f" % (T.sum(), T.sum() / n_in))

    P("=" * 70)
    t = T[0]
    try:
        tin = np.real(np.asarray(neff_in)).reshape(-1)
        tout = np.real(np.asarray(neff_out)).reshape(-1)
        P("input port modes  : %s" % ", ".join("m%d=%.4f" % (i, v) for i, v in enumerate(tin)))
        P("output port modes : %s" % ", ".join("m%d=%.4f" % (i, v) for i, v in enumerate(tout)))
        P("expect @x=105 (local FDE/FDTD): TE0~2.71, TE1~2.25")
        cy_out = ycenter(md["modes_out"], "Ey")
        i_te1 = int(np.argmin(np.abs(tin - 2.248)))
        o_bot = int(np.argmin(np.abs(np.asarray(cy_out) - (-0.875)))) if cy_out is not None else None
        o_top = int(np.argmin(np.abs(np.asarray(cy_out) - (-0.10)))) if cy_out is not None else None
        P("IDENT BY NEFF/CENTROID: i_TE1=%d (neff=%.4f) | o_bottom=%s | o_top=%s"
          % (i_te1, tin[i_te1], o_bot, o_top))
        if o_bot is not None and o_top is not None:
            P("KEY  in%d (TE1) -> out%d (bottom TE0) = %.5f   <== device function"
              % (i_te1, o_bot, t[o_bot, i_te1]))
            P("     in%d (TE1) -> out%d (top TE0)    = %.5f   (unwanted)"
              % (i_te1, o_top, t[o_top, i_te1]))
            P("     column sum of in%d = %.5f  (<=1 ; ==1.0000 means no radiation loss inside the basis)"
              % (i_te1, t[:, i_te1].sum()))
    except Exception as e:                                      # noqa: BLE001
        P("ident failed: %r" % (e,))

    flush()
    P("[done] log written to %s" % LOG)
    flush()
except Exception:                                              # noqa: BLE001
    P("TRACEBACK:")
    P(traceback.format_exc())
    flush()
    raise
