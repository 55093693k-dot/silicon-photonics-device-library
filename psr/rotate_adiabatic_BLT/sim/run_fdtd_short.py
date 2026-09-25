# run_fdtd_short.py -- 公开示例的做法（3D FDTD + ModeSource + ModeMonitor）的【缩短版】决定性实验
# 目的：只验证 bi-level taper 段（x = 0…105 µm）的 **TM0 -> TE1** 转换（该示例机理的表征点）
#       若该转换成立 ⇒ 器件机理在本平台成立，EME 的坏结果 = 方法（有限模式基）问题
# 用法：python run_fdtd_short.py            -> dry-run（只报成本，0 FlexCredit）
#       python run_fdtd_short.py --submit   -> 提交（浅 batch：TE + TM 两个源）
import os
import sys

import numpy as np
import tidy3d as td
import tidy3d.web as web

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_bilevel_psr as rbp

P = rbp.P
LAM = 1.55
X_IN = -4.0                     # 输入模式源位置（位于 0.45 µm 输入直波导内）
X_MON = P["L_blt"] + P["L_s"]   # = 105 µm：该示例机理的表征点（此处 mode0=TE0, mode1=TE1）
X0, X1 = -6.0, 109.0            # 仿真域 x 范围（含 taper 100 + 直段 5 + 余量）
Y_SPAN, Z_SPAN = 3.2, 2.2       # y 3.2：taper 横向 ±0.83 外留 >λ/2 的 PML 余量（易错项）

freq = td.C_0 / LAM
# z 向保留已验证的 20 nm 分辨率（避免"粗网格出结论"）；面内用 37 nm 降本
_ov = td.MeshOverrideStructure(
    geometry=td.Box(center=((X0 + X1) / 2, 0.0, 0.11), size=(X1 - X0, Y_SPAN, 1.0)),
    dl=(0.037, 0.037, 0.02),
)
common = dict(
    grid_spec=td.GridSpec.auto(min_steps_per_wvl=12, wavelength=LAM, override_structures=(_ov,)),
    structures=rbp.make_structures(P),
    boundary_spec=td.BoundarySpec.all_sides(boundary=td.PML()),
    medium=td.material_library["SiO2"]["Palik_Lossless"],
    run_time=3.6e-12,          # 约束：L=109 µm, ng≈4.5 -> t_transit≈1.6 ps，run_time ≥ 2×
)
mon = td.ModeMonitor(center=(X_MON, 0.0, 0.11), size=(0.0, Y_SPAN, Z_SPAN),
                     freqs=[freq], mode_spec=td.ModeSpec(num_modes=2), name="taper_out")


def build(idx):
    src = td.ModeSource(
        center=(X_IN, 0.0, 0.11), size=(0.0, Y_SPAN, Z_SPAN), direction="+",
        source_time=td.GaussianPulse(freq0=freq, fwidth=freq / 12),
        mode_spec=td.ModeSpec(num_modes=2), mode_index=idx,
    )
    return td.Simulation(center=((X0 + X1) / 2, 0.0, 0.0), size=(X1 - X0, Y_SPAN, Z_SPAN),
                         sources=[src], monitors=[mon], **common)


s0 = build(0)     # TE0 注入
s0.validate_pre_upload()
s1 = build(1)     # TM0 注入
s1.validate_pre_upload()
print("=== shortened FDTD: bi-level taper only (x = %.0f..%.0f um) ===" % (X0, X1))
print("mode source x=%.1f | mode monitor x=%.1f | grid=20 steps/lambda | all-side PML" % (X_IN, X_MON))
print("mesh cells (approx):", int(np.prod([c / 0.02 for c in s1.grid_spec.grid_size(s1)]) if False else 0) or "n/a")

job_te = web.Job(simulation=s0, task_name="psr_taper_short_TE", verbose=False)
job_tm = web.Job(simulation=s1, task_name="psr_taper_short_TM", verbose=False)
POL = sys.argv[sys.argv.index("--pol") + 1] if "--pol" in sys.argv else "tm"   # tm | te | both
jobs = {}
if POL in ("te", "both"):
    jobs["TE"] = job_te
if POL in ("tm", "both"):
    jobs["TM"] = job_tm
print("polarization:", POL, "| jobs:", list(jobs))
tot = 0.0
for k, j in jobs.items():
    c = web.estimate_cost(j.task_id)
    tot += c
    print("estimate %s (FlexCredit): %s" % (k, c))
print("estimate TOTAL (FlexCredit): %s" % tot)

if "--submit" not in sys.argv:
    print("[dry-run] not submitted")
    sys.exit(0)
for k, j in jobs.items():
    res = j.run(path=os.path.join(HERE, "data_taper_%s.hdf5" % k))
    amp = res["taper_out"].amps.sel(direction="+")
    T = (np.abs(np.asarray(amp.values).squeeze()) ** 2)
    print("%s input  -> T(mode0)=%.4f  T(mode1)=%.4f" % (k, T.ravel()[0], T.ravel()[1]))
