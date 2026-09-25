# run_eme_passive.py -- 最小付费变体：EME(3D) constraint="passive" + 单波长 1.55 um
# 目的：检查含辐射/回反射损耗的**绝对插损**能否由 EME 给出
#       （`unitary` 会强制 S 矩阵酉性，把基内列和压到 1，看不到真实损耗通道）
# 用法：python run_eme_passive.py [--modes 2]  -> dry-run（0 FlexCredit，只报成本）
#       python run_eme_passive.py [--modes 2] --submit   -> 真正提交
import os
import sys

import tidy3d as td
from tidy3d import web

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import run_bilevel_psr as rbp                     # 复用已验证的几何/EME 配方

P = rbp.P
LAM = 1.55                                        # um
nmodes = int(sys.argv[sys.argv.index("--modes") + 1]) if "--modes" in sys.argv else 2
sim = rbp.make_eme_sim(P, num_modes=nmodes)
sim = sim.updated_copy(constraint="passive", freqs=[td.C_0 / LAM])
sim.validate_pre_upload()
print("=== variant: passive + num_modes=%d + single wavelength %.2f um ===" % (nmodes, LAM))
print("constraint =", sim.constraint, "| n_freqs =", len(sim.freqs),
      "| eme cells =", sim.eme_grid_spec.boundaries.size + 1)

task = "bilevel_psr_eme_passive_m%d_lam1550" % nmodes
job = web.Job(simulation=sim, task_name=task, verbose=False)
print("estimate_cost (FlexCredit):", web.estimate_cost(job.task_id))
if "--submit" not in sys.argv:
    print("[dry-run] not submitted")
    sys.exit(0)
data = job.run(path=os.path.join(HERE, "data_%s.hdf5" % task))
S = data.smatrix.S21
print("S21 dims:", tuple(S.dims), "shape:", S.shape)
print("data file:", os.path.join(HERE, "data_%s.hdf5" % task))
