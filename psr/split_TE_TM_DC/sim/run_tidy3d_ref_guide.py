# -*- coding: utf-8 -*-
"""
run_tidy3d_ref_guide.py —— 直波导参照测试（隔离 TM 注入/监测问题）
几何: 单根 450x220nm SOI 条（无耦合臂），x 向 10um + 2um stub
运行: python run_tidy3d_ref_guide.py [--mode BOTH]
期望: TE (mode0) 与 TM (mode1) 传输都应 ~1.0（|amp|^2），否则源/监测设置有误
"""
from __future__ import annotations
import argparse, os
import numpy as np
import tidy3d as td
import tidy3d.web as web

OUT = os.path.dirname(os.path.abspath(__file__))
N_SI, N_SIO2 = 3.476, 1.444
H, W = 0.22, 0.45
L = 10.0
STUB = 2.0
X0, X1 = -STUB, L + STUB
FREQ0 = 299792458.0 / 1.55e-6
FREQW = 6.0e12
RUN_TIME = 2.0e-12
MESH = (0.02, 0.01, 0.01)


def box(cx, cy, cz, sx, sy, sz):
    return td.Box(center=(cx, cy, cz), size=(sx, sy, sz))


def build(mode_index):
    si = td.Medium(permittivity=N_SI ** 2)
    sio2 = td.Medium(permittivity=N_SIO2 ** 2)
    simx = X1 - X0
    cx = (X0 + X1) / 2
    structs = [
        td.Structure(geometry=box(cx, 0, 0, simx, 3.0, 3.0), medium=sio2),
        td.Structure(geometry=box(cx, H / 2, 0, simx + 2.4, H, W), medium=si),
    ]
    mo = td.ModeSpec(num_modes=2, target_neff=2.0)
    freqs = np.linspace(FREQ0 - FREQW, FREQ0 + FREQW, 7)
    src = [td.ModeSource(center=(X0 + 0.3, H / 2, 0), size=(0, 1.0, 1.0),
                         source_time=td.GaussianPulse(freq0=FREQ0, fwidth=FREQW),
                         mode_index=mode_index, mode_spec=mo, direction="+")]
    mon = td.ModeMonitor(center=(X1 - 0.5, H / 2, 0), size=(0, 1.0, 1.0), freqs=freqs,
                         mode_spec=mo, name="out", store_fields_direction="+")
    ov = td.MeshOverrideStructure(geometry=box(cx, 0.11, 0, simx, 0.5, 1.4), dl=MESH)
    return td.Simulation(center=(cx, 0, 0), size=(simx, 2.4, 4.0),
                         grid_spec=td.GridSpec(override_structures=[ov], wavelength=1.55),
                         structures=structs, sources=src, monitors=[mon],
                         run_time=RUN_TIME, boundary_spec=td.BoundarySpec.all_sides(td.PML()))


def run(mode_index, ptag):
    sim = build(mode_index)
    name = "REF_GUIDE_%s" % ptag
    sd = web.run(sim, task_name=name, path=os.path.join(OUT, name + ".hdf5"), verbose=False)
    out = sd["out"]
    amps = np.abs(out.amps.sel(direction="+").values)
    neff = np.abs(out.n_eff.values)
    lam = 299792458.0 / np.array(out.amps.f) * 1e9
    with open(os.path.join(OUT, "tidy3d_REF_%s.csv" % ptag), "w") as f:
        f.write("lambda_nm,mode0_amp,mode1_amp,neff0,neff1\n")
        for k in range(len(lam)):
            f.write("%.2f,%.6e,%.6e,%.4f,%.4f\n" %
                    (lam[k], amps[k, 0], amps[k, 1],
                     neff[k, 0] if neff.shape[1] > 0 else 0,
                     neff[k, 1] if neff.shape[1] > 1 else 0))
    print(ptag, "amps:", amps.max(axis=0).round(4), "neff@1550:", neff[int(np.argmin(np.abs(lam-1550)))].round(4))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["TE", "TM", "BOTH"], default="BOTH")
    a = ap.parse_args()
    if a.mode in ("TE", "BOTH"):
        run(0, "TE")
    if a.mode in ("TM", "BOTH"):
        run(1, "TM")
