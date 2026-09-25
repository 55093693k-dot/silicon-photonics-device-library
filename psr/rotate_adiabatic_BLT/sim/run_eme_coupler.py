# run_eme_coupler.py -- 耦合器段（绝热段）EME 验证
# 域 x 100…410 µm；端口严格放在该示例的耦合器边界 x=105（入口，TE1 注入）/ x=405（出口）
# 单频 λ=1.55 µm + constraint="passive"（要绝对 IL，不能 unitary）+ 存端口模式（按 neff 认模）
# 用法：python run_eme_coupler.py [--modes 4] [--y 6] [--submit]   （不带 --submit = dry-run，0 FlexCredit）
import os
import sys

import numpy as np
import tidy3d as td
from tidy3d import web

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_bilevel_psr as rbp                              # 复用已验证几何/EME 配方

P = rbp.P
LAM = 1.55
NMODES = int(sys.argv[sys.argv.index("--modes") + 1]) if "--modes" in sys.argv else 4
YS = float(sys.argv[sys.argv.index("--y") + 1]) if "--y" in sys.argv else 6.0

XP1 = P["L_blt"] + P["L_s"]                                # 105 µm：junction / 耦合器入口（该示例：mode1=TE1）
XP2 = XP1 + P["L_ac"]                                      # 405 µm：耦合器出口（该示例：mode3→底部 TE0）
X0, X1 = XP1 - 5.0, XP2 + 5.0                              # 域 100…410 µm
PAD = 5.0

# EME 网格（照抄已验证分区法）：junction 附近 0.5 µm/格 + 耦合器 5 µm/格（绝热）
# 注意（源码 eme/simulation.py:807）：EME 网格区间 = 域**扣除 port_offsets**，
#   端口面（105/405）自动成为首/末边界 ⇒ 用户 boundaries 必须**严格落在两端口之间**（易错项）。
xs = np.concatenate([np.arange(XP1 + 0.5, XP1 + PAD + 1e-9, 0.5),
                     np.arange(XP1 + PAD, XP2 + 1e-9, 5.0)])
bnd = sorted({round(float(v), 5) for v in xs if XP1 + 1e-6 < float(v) < XP2 - 1e-6})
grid = td.EMEExplicitGrid(mode_specs=[td.EMEModeSpec(num_modes=NMODES) for _ in range(len(bnd) + 1)],
                          boundaries=tuple(bnd))
assert all(abs(v - XP1) > 1e-6 and abs(v - XP2) > 1e-6 for v in bnd), \
    "端口面由 port_offsets 隐式给定，不得作为显式 boundary（否则 ValidationError）"

sim0 = rbp.make_eme_sim(P, num_modes=NMODES)               # 已验证配方 → 只改域/网格/约束/频率/监视器
sim = sim0.updated_copy(
    center=((X0 + X1) / 2, 0.0, 0.11), size=(X1 - X0, YS, 1.75),
    eme_grid_spec=grid, constraint="passive", freqs=[td.C_0 / LAM],
    port_offsets=(PAD, PAD),
    monitors=[td.EMEModeSolverMonitor(name="modes_in", num_modes=NMODES,
                                      size=(0, td.inf, td.inf), center=(XP1, 0.0, 0.11)),
              td.EMEModeSolverMonitor(name="modes_out", num_modes=NMODES,
                                      size=(0, td.inf, td.inf), center=(XP2, 0.0, 0.11)),
              td.EMECoefficientMonitor(name="coeffs", size=(td.inf, td.inf, td.inf))],
)
sim.validate_pre_upload()
print("=== coupler EME: 域 x %.0f…%.0f µm | 端口 x=%.0f(入口) / %.0f(出口) ===" % (X0, X1, XP1, XP2))
print("modes/port=%d | EME cells=%d (%d internal bnd; 0.5µm@junction + 5µm@coupler) | y=%.1f | f=%d | constraint=%s"
      % (NMODES, len(bnd) + 1, len(bnd), YS, len(sim.freqs), sim.constraint))

task = "psr_coupler_eme_m%d_y%g" % (NMODES, YS)
job = web.Job(simulation=sim, task_name=task, verbose=False)
print("estimate_cost (FlexCredit):", web.estimate_cost(job.task_id))
if "--submit" not in sys.argv:
    print("[dry-run] not submitted")
    sys.exit(0)

data = job.run(path=os.path.join(HERE, "data_%s.hdf5" % task))
S = data.smatrix.S21
print("S21 dims:", tuple(S.dims), "shape:", S.shape)
for nm, tag in (("modes_in", "入口 x=105"), ("modes_out", "出口 x=405")):
    try:
        md = getattr(data, nm)
        print("neff %s (%s):" % (nm, tag), np.round(np.real(np.atleast_2d(md.n_complex.values)[0]), 4))
    except Exception as e:                                  # noqa: BLE001
        print("neff %s: 读取失败 %s" % (nm, e))
print("|S21|^2 (out x in) =")
print(np.round(np.abs(S.values) ** 2, 4))
print("data file:", os.path.join(HERE, "data_%s.hdf5" % task))
