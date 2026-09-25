# -*- coding: utf-8 -*-
"""
run_tidy3d_splitter.py —— TE/TM 模式分离器 3D FDTD 验证（Tidy3D）
几何(FDE 最优): w1=w2=450nm strip, gap=150nm, L=31.5um (+2um access)
运行: python run_tidy3d_splitter.py [--source TE|TM|BOTH]
判读: out_wg1/out_wg2 各 mode 幅值 -> TE 留 wg1, TM 交叉到 wg2
备注: 该方案为等宽 DC + 选择长度(偏振相关耦合长度比≈3.5)；需 3D 验证 Δn 标定。
"""
from __future__ import annotations
import argparse, os
import numpy as np
import tidy3d as td
import tidy3d.web as web

OUT = os.path.dirname(os.path.abspath(__file__))
N_SI, N_SIO2 = 3.476, 1.444
H = 0.22
W1 = W2 = 0.45
GAP = 0.15
L = 31.5
STUB = 2.0
X0, X1 = -STUB, L + STUB
FREQ0 = 299792458.0 / 1.55e-6
FREQW = 6.0e12
RUN_TIME = 3.0e-12
MESH = (0.02, 0.01, 0.01)   # 10nm: mesh-convergence verified


def box(cx, cy, cz, sx, sy, sz):
    return td.Box(center=(cx, cy, cz), size=(sx, sy, sz))


def build(source):
    si = td.Medium(permittivity=N_SI ** 2)
    sio2 = td.Medium(permittivity=N_SIO2 ** 2)
    simx = X1 - X0
    cx = (X0 + X1) / 2
    z1 = -GAP / 2 - W1 / 2
    z2 = GAP / 2 + W2 / 2
    # wg2: S-bend access (offset 1.35um at x=0 -> coupling position at RAMP)
    OFF = 1.35
    RAMP = 3.0
    verts = [(X1 + 1.2, z2 + W2 / 2), (X1 + 1.2, z2 - W2 / 2), (RAMP, z2 - W2 / 2),
             (0.0, z2 + OFF - W2 / 2), (0.0, z2 + OFF + W2 / 2), (RAMP, z2 + W2 / 2)]
    structs = [
        td.Structure(geometry=box(cx, 0, 0, simx, 4.8, 6.0), medium=sio2),
        td.Structure(geometry=box(cx, H / 2, z1, simx + 2.4, H, W1), medium=si),
        td.Structure(geometry=td.PolySlab(vertices=verts, axis=1, slab_bounds=(0.0, H)), medium=si),
    ]
    mo = td.ModeSpec(num_modes=2, target_neff=2.0)
    freqs = np.linspace(FREQ0 - FREQW, FREQ0 + FREQW, 7)
    src = []
    if source in ("TE", "BOTH"):
        src.append(td.ModeSource(center=(X0 + 0.3, H / 2, z1), size=(0, 1.0, 1.0),
                                 source_time=td.GaussianPulse(freq0=FREQ0, fwidth=FREQW),
                                 mode_index=0, mode_spec=mo, direction="+"))
    if source in ("TM", "BOTH"):
        src.append(td.ModeSource(center=(X0 + 0.3, H / 2, z1), size=(0, 1.0, 1.0),
                                 source_time=td.GaussianPulse(freq0=FREQ0, fwidth=FREQW),
                                 mode_index=1, mode_spec=mo, direction="+"))
    mons = [
        td.ModeMonitor(center=(X1 - 0.4, H / 2, z1), size=(0, 1.0, 0.60), freqs=freqs,
                       mode_spec=mo, name="out_wg1", store_fields_direction="+"),
        td.ModeMonitor(center=(X1 - 0.4, H / 2, z2), size=(0, 1.0, 0.60), freqs=freqs,
                       mode_spec=mo, name="out_wg2", store_fields_direction="+"),
    ]
    ov = td.MeshOverrideStructure(geometry=box(cx, 0.11, 0, simx, 0.5, 2.0), dl=MESH)
    return td.Simulation(center=(cx, 0, 0), size=(simx, 2.8, 5.0),
                         grid_spec=td.GridSpec(override_structures=[ov], wavelength=1.55),
                         structures=structs, sources=src, monitors=mons,
                         run_time=RUN_TIME, boundary_spec=td.BoundarySpec.all_sides(td.PML()))


def run(source):
    sim = build(source)
    name = "SPLITTER_te_tm_%s" % source
    sd = web.run(sim, task_name=name, path=os.path.join(OUT, name + ".hdf5"), verbose=True)
    for mon in ("out_wg1", "out_wg2"):
        out = sd[mon]
        amps = np.abs(out.amps.sel(direction="+").values)
        lam = 299792458.0 / np.array(out.amps.f) * 1e9
        with open(os.path.join(OUT, "tidy3d_SPLIT_%s_%s.csv" % (source, mon)), "w") as f:
            f.write("lambda_nm," + ",".join("mode%d" % i for i in range(amps.shape[1])) + "\n")
            for l, row in zip(lam, amps):
                f.write("%.2f," % l + ",".join("%.6e" % v for v in row) + "\n")
        print(mon, amps.max(axis=0).round(4))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", choices=["TE", "TM", "BOTH"], default="BOTH")
    a = ap.parse_args()
    for s in (["TE", "TM"] if a.source == "BOTH" else [a.source]):
        run(s)
    print("done: TE->wg1 mode0 | TM->wg2 mode0 expected")
