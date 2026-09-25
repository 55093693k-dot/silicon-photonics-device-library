# -*- coding: utf-8 -*-
"""plot_port_mode_fields.py -- 离线画「耦合器段端口模式真实场分布」（0 FlexCredit）

输入：run_eme_coupler.py 已下载的 hdf5（默认取项目目录里的 6 模式跑）
输出：../data/PORT_MODE_FIELDS.png  —— 入口(x=105um)/出口(x=405um) 各模式 |Ey|^2 面图 + neff

坐标（本 EME 设置）：传播 = x；横向（两臂分开方向）= y；垂直 = z。
本脚本按**维度名**选取（不假设轴序），自动裁到模式有场的横向范围以利观察。

用法:  python plot_port_mode_fields.py [hdf5 路径]
依赖:  tidy3d, numpy, matplotlib（本地免费；不提交任何任务）
"""
import os
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import tidy3d as td  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_H5 = os.path.join(HERE, "..", "..", "..", "..", "PSR", "bilevel_psr",
                          "data_psr_coupler_eme_m6_y6.hdf5")
OUT = os.path.join(HERE, "..", "data", "PORT_MODE_FIELDS.png")


def load(path):
    for meth in ("from_hdf5", "from_file"):
        f = getattr(td.EMESimulationData, meth, None)
        if f is None:
            continue
        try:
            d = f(path)
            print("loaded via EMESimulationData.%s()" % meth)
            return d
        except Exception as e:  # noqa: BLE001
            print("%s failed: %r" % (meth, str(e)[:200]))
    return None


def mode_plane(ds, comp="Ey"):
    """返回 (w[mode, z, y], y 轴, z 轴, neff)。按维度名选取，去掉所有长度 1 的轴。"""
    c = getattr(ds, comp, None)
    if c is None:
        return None
    for d in list(c.dims):
        if c.sizes[d] == 1 and d != "mode_index":
            c = c.isel({d: 0})
    w = np.abs(np.asarray(c.values)) ** 2
    dims = list(c.dims)
    if w.ndim == 2:                       # 单模式：补一个 mode 轴
        w = w[None, :, :]
        dims = ["mode_index"] + dims
    order = [dims.index(a) for a in ("mode_index", "z", "y")]
    w = np.transpose(w, order)            # -> (mode, z(vertical), y(lateral))
    neff = np.real(np.asarray(ds.n_complex.values)).reshape(-1)
    return w, np.asarray(c.coords["y"].values), np.asarray(c.coords["z"].values), neff


def main():
    path = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else DEFAULT_H5)
    print("hdf5 exists: %s -> %s" % (os.path.exists(path), path))
    if not os.path.exists(path):
        sys.exit(1)
    data = load(path)
    if data is None:
        sys.exit(1)
    md = data.monitor_data
    print("monitor_data keys: %s" % list(md.keys()))

    keys = [k for k in ("modes_in", "modes_out") if k in md]
    panels = [(k, mode_plane(md[k])) for k in keys]
    panels = [(k, p) for k, p in panels if p is not None]
    if not panels:
        print("no usable panels")
        sys.exit(1)
    for k, p in panels:
        print("[%s] |Ey|^2 shape=%s neff=%s" % (k, p[0].shape, np.round(p[3], 4)))

    # 横向裁剪：所有模式合计有场（>1% 峰值）的 y 范围 + 20% 余量
    y = panels[0][1][1]
    prof = np.zeros_like(y)
    for _, (w, yy, zz, neff) in panels:
        prof = prof + w.sum(axis=(0, 1))
    keep = np.where(prof > 0.01 * prof.max())[0]
    margin = max(2, int(0.20 * (keep[-1] - keep[0] + 1)))
    i0 = max(0, keep[0] - margin)
    i1 = min(len(y), keep[-1] + 1 + margin)
    print("lateral crop: y = %.3f .. %.3f um (full %.3f .. %.3f)"
          % (y[i0], y[i1 - 1], y[0], y[-1]))

    nmax = max(p[0].shape[0] for _, p in panels)
    fig, axs = plt.subplots(len(panels), nmax,
                            figsize=(2.5 * nmax + 0.8, 3.0 * len(panels) + 1.0), squeeze=False)
    for r, (key, (w, yy, zz, neff)) in enumerate(panels):
        for c in range(nmax):
            ax = axs[r][c]
            if c < w.shape[0]:
                im = ax.imshow(w[c][:, i0:i1], origin="lower", aspect="auto", cmap="inferno",
                               extent=[float(yy[i0]), float(yy[i1 - 1]), float(zz[0]), float(zz[-1])])
                fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
                ax.set_title("%s m%d  neff=%.4f" % (key, c, neff[c]), fontsize=8)
                ax.set_xlabel("lateral y (um)", fontsize=7)
            else:
                ax.axis("off")
            if c == 0:
                ax.set_ylabel("vertical z (um)", fontsize=7)
            ax.tick_params(labelsize=6)
    fig.suptitle("Port mode fields |Ey|^2 @1.55um -- adiabatic coupler section\n"
                 "(in: x=105um, out: x=405um; 6 modes, passive)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT, dpi=130)
    print("saved: %s (%d B)" % (os.path.abspath(OUT), os.path.getsize(OUT)))


if __name__ == "__main__":
    main()
