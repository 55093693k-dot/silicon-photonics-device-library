# -*- coding: utf-8 -*-
"""
run_tidy3d_splitter_field.py —— 用场积分法测 TE/TM 分离（不依赖端口模式序号）
输出面 FieldMonitor(x-normal)，按臂的 z 窗口积分 x 向能流 Sx=0.5*Re(Ey*Hz*-Ez*Hy*)
并给出每臂 Ex/Ey 占比（判别偏振）；输入接入段仅 wg1（干净注入 TE0/TM0）
运行: python run_tidy3d_splitter_field.py [--source TE|TM|BOTH]
"""
from __future__ import annotations
import argparse, os
import numpy as np
import tidy3d as td
import tidy3d.web as web

OUT = os.path.dirname(os.path.abspath(__file__))
N_SI, N_SIO2 = 3.476, 1.444
H, W1, W2, GAP, L = 0.22, 0.45, 0.45, 0.15, 31.5
STUB = 2.0
RAMP = 3.0
X0, X1 = -STUB, L + STUB
FREQ0 = 299792458.0 / 1.55e-6
FREQW = 6.0e12
RUN_TIME = 3.0e-12
MESH = (0.02, 0.01, 0.01)
Z1 = -GAP / 2 - W1 / 2
Z2 = GAP / 2 + W2 / 2


def box(cx, cy, cz, sx, sy, sz):
    return td.Box(center=(cx, cy, cz), size=(sx, sy, sz))


def build(src_mode):
    si = td.Medium(permittivity=N_SI ** 2)
    sio2 = td.Medium(permittivity=N_SIO2 ** 2)
    simx = X1 - X0
    cx = (X0 + X1) / 2
    # wg2: S-bend access (away by 1.35um at x=0, ramps to coupling position over RAMP)
    OFF = 1.35
    verts = [(X1 + 1.2, Z2 + W2 / 2), (X1 + 1.2, Z2 - W2 / 2), (RAMP, Z2 - W2 / 2),
             (0.0, Z2 + OFF - W2 / 2), (0.0, Z2 + OFF + W2 / 2), (RAMP, Z2 + W2 / 2)]
    wg2 = td.PolySlab(vertices=verts, axis=1, slab_bounds=(0.0, H))
    structs = [
        td.Structure(geometry=box(cx, 0, 0, simx, 4.8, 6.0), medium=sio2),
        td.Structure(geometry=box(cx, H / 2, Z1, simx + 2.4, H, W1), medium=si),
        td.Structure(geometry=wg2, medium=si),
    ]
    mo = td.ModeSpec(num_modes=2, target_neff=2.0)
    src = [td.ModeSource(center=(X0 + 0.3, H / 2, Z1), size=(0, 1.0, 1.0),
                         source_time=td.GaussianPulse(freq0=FREQ0, fwidth=FREQW),
                         mode_index=src_mode, mode_spec=mo, direction="+")]
    fm = td.FieldMonitor(center=(X1 - 0.4, 0.11, 0.0), size=(0.0, 1.0, 3.0),
                         freqs=[FREQ0], name="out_plane")
    fmi = td.FieldMonitor(center=(X0 + 0.5, 0.11, 0.0), size=(0.0, 1.0, 2.0),
                          freqs=[FREQ0], name="in_plane")
    ov = td.MeshOverrideStructure(geometry=box(cx, 0.11, 0, simx, 0.5, 2.0), dl=MESH)
    return td.Simulation(center=(cx, 0, 0), size=(simx, 2.8, 6.0),
                         grid_spec=td.GridSpec(override_structures=[ov], wavelength=1.55),
                         structures=structs, sources=src, monitors=[fm, fmi],
                         run_time=RUN_TIME, boundary_spec=td.BoundarySpec.all_sides(td.PML()))


def analyze(sd, tag):
    fd = sd["out_plane"]
    Sx = 0.5 * np.real(fd.Ey * np.conj(fd.Hz) - fd.Ez * np.conj(fd.Hy))
    pEx = np.abs(fd.Ex) ** 2
    pEy = np.abs(fd.Ey) ** 2
    lines = ["port,power_Sx,Ex_frac,Ey_frac",
             "dims,%s,," % (",".join(Sx.dims))]
    tot = 0.0
    _pf = {}
    for key, zc in (("wg1", Z1), ("wg2", Z2)):
        sl = slice(zc - 0.45, zc + 0.45)
        P = float(Sx.sel(z=sl).sum().values)
        aE = float(pEx.sel(z=sl).sum().values)
        aY = float(pEy.sel(z=sl).sum().values)
        d = aE + aY + 1e-30
        lines.append("%s,%.6e,%.4f,%.4f" % (key, P, aE / d, aY / d))
        tot += P
        _pf[key] = P
    lines.append("wg1+wg2,%.6e,," % tot)
    lines.append("frac_wg1,%.4f,," % (_pf["wg1"] / tot))
    lines.append("frac_wg2,%.4f,," % (_pf["wg2"] / tot))
    if "in_plane" in sd:
        fi = sd["in_plane"]
        Sxi = 0.5 * np.real(fi.Ey * np.conj(fi.Hz) - fi.Ez * np.conj(fi.Hy))
        Pin = float(Sxi.sel(z=slice(-1.0, 1.0)).sum().values)
        lines.append("P_in,%.6e,," % Pin)
        for key, zc in (("wg1", Z1), ("wg2", Z2)):
            sl = slice(zc - 0.45, zc + 0.45)
            Ti = float(Sx.sel(z=sl).sum().values) / Pin
            lines.append("T_%s,%.4f,," % (key, Ti))
        lines.append("T_total,%.4f,," % (tot / Pin))
    with open(os.path.join(OUT, "tidy3d_FIELD_%s.csv" % tag), "w") as f:
        f.write("\n".join(lines))
    print("\n".join(lines))


def run(src_mode, tag):
    sim = build(src_mode)
    name = "SPLIT_FIELD_%s" % tag
    sd = web.run(sim, task_name=name, path=os.path.join(OUT, name + ".hdf5"), verbose=False)
    analyze(sd, tag)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", choices=["TE", "TM", "BOTH"], default="BOTH")
    ap.add_argument("--hdf5", default=None, help="仅分析已下载的 hdf5（不重新上云）")
    ap.add_argument("--tag", default="TE", help="配合 --hdf5 使用的输出标签")
    a = ap.parse_args()
    if a.hdf5:
        analyze(td.SimulationData.from_file(a.hdf5), a.tag)
    else:
        if a.source in ("TE", "BOTH"):
            run(0, "TE")
        if a.source in ("TM", "BOTH"):
            run(1, "TM")
