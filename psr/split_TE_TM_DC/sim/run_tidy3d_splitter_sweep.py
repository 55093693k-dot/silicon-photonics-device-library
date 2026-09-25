# -*- coding: utf-8 -*-
"""
run_tidy3d_splitter_sweep.py —— TE/TM 分离器 Tidy3D 内重定标批量扫描
配置: gap in {200,300,400,500}nm (L=31.5um)  +  L in {10,20,40,60}um (gap=150nm)
每配置跑 TE 与 TM 两次 -> 判定 TE 直通(wg1)>95% 且 TM 交叉(wg2)>95%
运行: python run_tidy3d_splitter_sweep.py [--which gap|L|both] [--mode TE|TM|BOTH]
输出: sweep_results.csv / sweep_summary.md / tidy3d_SWEEP_<cfg>_<pol>.csv
注意: 16 次云任务(默认 both+BOTH)，可先用 --mode TM 只跑关键通道省额度
"""
from __future__ import annotations
import argparse, os, csv
import numpy as np
import tidy3d as td
import tidy3d.web as web

OUT = os.path.dirname(os.path.abspath(__file__))
N_SI, N_SIO2 = 3.476, 1.444
H = 0.22
W1 = W2 = 0.45
STUB = 2.0
FREQ0 = 299792458.0 / 1.55e-6
FREQW = 6.0e12
RUN_TIME = 3.0e-12
MESH = (0.02, 0.01, 0.01)     # 10nm (mesh-convergence verified)
NUM_MODES = 2
NEFF_TE0, NEFF_TM0 = 2.30, 1.72   # 用于按 neff 识别 TE0 / TM0


def box(cx, cy, cz, sx, sy, sz):
    return td.Box(center=(cx, cy, cz), size=(sx, sy, sz))


def build(gap, L, source):
    si = td.Medium(permittivity=N_SI ** 2)
    sio2 = td.Medium(permittivity=N_SIO2 ** 2)
    X0, X1 = -STUB, L + STUB
    simx = X1 - X0
    cx = (X0 + X1) / 2
    z1 = -gap / 2 - W1 / 2
    z2 = gap / 2 + W2 / 2
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
    mo = td.ModeSpec(num_modes=NUM_MODES, target_neff=2.0)
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


def _pow_by_neff(mon_data, target_neff):
    """power |amp|^2 over freq for the mode whose neff ~ target (label-safe)."""
    amps = mon_data.amps.sel(direction="+").values
    neff = np.abs(mon_data.n_eff.values)
    p = np.abs(amps) ** 2
    idx = np.argmin(np.abs(neff - target_neff), axis=1)
    return np.array([p[i, idx[i]] for i in range(p.shape[0])])


def run_cfg(gap_nm, L_um, mode):
    gap, L = gap_nm * 1e-3, float(L_um)
    tag = "g%03d_L%04.1f" % (gap_nm, L)
    sim = build(gap, L, mode)
    name = "SPLIT_SWEEP_%s_%s" % (tag, mode)
    sd = web.run(sim, task_name=name, path=os.path.join(OUT, name + ".hdf5"), verbose=True)
    lam = 299792458.0 / np.array(sd["out_wg1"].amps.f) * 1e9
    res = {"tag": tag, "gap_nm": gap_nm, "L_um": L, "pol": mode}
    for mon, key in (("out_wg1", "wg1"), ("out_wg2", "wg2")):
        pTE = _pow_by_neff(sd[mon], NEFF_TE0)
        pTM = _pow_by_neff(sd[mon], NEFF_TM0)
        res[key + "_TE"] = pTE
        res[key + "_TM"] = pTM
        with open(os.path.join(OUT, "tidy3d_SWEEP_%s_%s_%s.csv" % (tag, mode, mon)), "w") as f:
            f.write("lambda_nm,TE0_pow,TM0_pow\n")
            for l, a, b in zip(lam, pTE, pTM):
                f.write("%.2f,%.6e,%.6e\n" % (l, a, b))
    i = int(np.argmin(np.abs(lam - 1550)))
    return res, i


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--which", choices=["gap", "L", "both"], default="both")
    ap.add_argument("--mode", choices=["TE", "TM", "BOTH"], default="BOTH")
    a = ap.parse_args()
    cfgs = []
    if a.which in ("gap", "both"):
        for g in (200, 300, 400, 500):
            cfgs.append((g, 31.5))
    if a.which in ("L", "both"):
        for L in (10, 20, 40, 60):
            cfgs.append((150, L))
    rows = []
    pols = ["TE", "TM"] if a.mode == "BOTH" else [a.mode]
    for (g, L) in cfgs:
        rec = {"gap_nm": g, "L_um": L}
        for pol in pols:
            r, i1550 = run_cfg(g, L, pol)
            if pol == "TE":
                rec["TE_bar_1550"] = float(r["wg1_TE"][i1550])
                rec["TE_bar_min"] = float(r["wg1_TE"].min())
                rec["TE_cross_1550"] = float(r["wg2_TE"][i1550])
            else:
                rec["TM_cross_1550"] = float(r["wg2_TM"][i1550])
                rec["TM_cross_min"] = float(r["wg2_TM"].min())
                rec["TM_bar_1550"] = float(r["wg1_TM"][i1550])
        ok = (rec.get("TE_bar_1550", 0) > 0.95) and (rec.get("TM_cross_1550", 0) > 0.95)
        rec["PASS"] = bool(ok)
        rows.append(rec)
        print("cfg g=%d L=%s -> TE_bar=%.3f TM_cross=%.3f PASS=%s"
              % (g, L, rec.get("TE_bar_1550", np.nan), rec.get("TM_cross_1550", np.nan), ok))
    keys = ["gap_nm", "L_um", "TE_bar_1550", "TE_bar_min", "TE_cross_1550",
            "TM_cross_1550", "TM_cross_min", "TM_bar_1550", "PASS"]
    with open(os.path.join(OUT, "sweep_results.csv"), "w", newline="") as f:
        w = csv.writer(f); w.writerow(keys)
        for r in rows:
            w.writerow([r.get(k, "") for k in keys])
    with open(os.path.join(OUT, "sweep_summary.md"), "w", encoding="utf-8") as f:
        f.write("# 分离器 Tidy3D 重定标汇总\n\n")
        f.write("| gap(nm) | L(um) | TE bar@1550 | TM cross@1550 | PASS |\n|---|---|---|---|---|\n")
        for r in rows:
            f.write("| %d | %s | %s | %s | %s |\n" % (
                r["gap_nm"], r["L_um"],
                ("%.3f" % r["TE_bar_1550"]) if "TE_bar_1550" in r else "-",
                ("%.3f" % r["TM_cross_1550"]) if "TM_cross_1550" in r else "-",
                "PASS" if r["PASS"] else "FAIL"))
        f.write("\n判定: TE bar(wg1)>95% 且 TM cross(wg2)>95%\n")
    print("saved sweep_results.csv / sweep_summary.md")


if __name__ == "__main__":
    main()
