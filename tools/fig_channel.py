# -*- coding: utf-8 -*-
"""fig_channel.py -- figure-channel labelling + PNG metadata reader.

WHY (the "figure channel rule")
-------------------------------
A figure's evidence grade depends on WHO drew it.  A Lumerical export is the
solver's own structure view / field / data plot; a matplotlib figure is a human
re-projection of the model.  The two must not be quoted as if equivalent, so:

  * official channels go first (index monitor / exportview);
  * matplotlib is FORBIDDEN unless (a) the official feature cannot make that
    figure and (b) the user approved it IN ADVANCE;
  * an approved external figure must say so NEXT TO THE FIGURE and in the
    report ("external plot, not an official export").

This module makes that requirement AUDITABLE instead of honour-system: it
renders the statement on the figure (a footer line) AND writes the same
statement into the PNG tEXt metadata, so a whole directory can be checked in
seconds without OCR (see fig_index.py).

Channel names (keep them identical everywhere):
  "official-monitor" : addindex() (material/index monitor) + image() + exportfigure()
  "official-view"    : exportview()
  "external"         : matplotlib / any non-Lumerical rendering  <- needs the label

Usage in a generator script
---------------------------
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import fig_channel as FC          # same directory on sys.path

    fig, ax = plt.subplots()
    ...draw...
    FC.save(fig, out_path, "external",
            "model top view, source sim/MODEL_OBJECTS.json")
    plt.close(fig)

`FC.save` = `fig.savefig` plus the footer label plus the tEXt chunk, so a script
only has to swap its `fig.savefig(...)` line (keep passing dpi/bbox_inches on).

No third-party dependency beyond matplotlib itself (the reader is a 20-line PNG
chunk walker, so verification works even if PIL is missing).
"""
from __future__ import annotations

import os
import struct

#: PNG tEXt key used by :func:`save` (and looked for by fig_index.py).
CHUNK = "External-Plot"

#: The ASCII statement per channel.  ASCII on purpose: PNG tEXt is Latin-1, and
#: keeping the metadata plain ASCII avoids any encoding surprises.
TEXT = {
    "official-monitor":
        "OFFICIAL export (Lumerical): addindex() + image() + exportfigure()",
    "official-view":
        "OFFICIAL export (Lumerical): exportview()",
    "external":
        "EXTERNAL PLOT (matplotlib) - NOT an official Lumerical export channel; "
        "prior user approval + this label required",
}

#: Chinese prefix for the ON-FIGURE footer of an external figure.
FIG_PREFIX = "外部绘图，非官方导出  /  "

#: Short ASCII form used ON THE FIGURE.  Kept short on purpose: a footer wider
#: than the canvas makes bbox_inches="tight" pad the saved image sideways (seen
#: once: a 1351 px figure came out 1561 px wide).  The full statement still
#: goes into the PNG metadata, which is what fig_index.py audits.
SHORT = "EXTERNAL PLOT - not an official export"

#: The phrase fig_index.py requires in an external figure's chunk (ASCII form of
#: the same statement) -- change it here and the tool follows.
MARK = "NOT an official Lumerical export channel"

#: Fonts tried for the on-figure footer, in order; first existing file wins.
#: (An explicit font FILE is used on purpose: matplotlib's *name*-based CJK
#: lookup is what is unreliable on this box, not the fonts.)
CJK_FONTS = (
    r"C:\Windows\Fonts\msyh.ttc",     # Microsoft YaHei
    r"C:\Windows\Fonts\simhei.ttf",   # SimHei
    r"C:\Windows\Fonts\simsun.ttc",   # SimSun
)


def _cjk_font(size=6.4):
    """FontProperties for the footer, or None when no CJK font file exists."""
    try:
        from matplotlib.font_manager import FontProperties
    except Exception:                                     # pragma: no cover
        return None
    for p in CJK_FONTS:
        if os.path.exists(p):
            try:
                return FontProperties(fname=p, size=size)
            except Exception:                             # pragma: no cover
                return None
    return None


def statement(channel="external", provenance=""):
    """The ASCII statement written to PNG metadata (and audited by fig_index)."""
    s = TEXT[channel]
    return s + (" | " + provenance if provenance else "")


def footer(channel="external", provenance=""):
    """The on-figure text: short Chinese+ASCII label (full text lives in metadata)."""
    if channel in TEXT and channel != "external":
        return TEXT[channel]
    return (FIG_PREFIX + SHORT) if _cjk_font() is not None else SHORT


def save(fig, path, channel="external", provenance="", **kw):
    """``fig.savefig(path, **kw)`` + the section 3.8 label (on-figure + metadata).

    ``kw`` is passed straight through, so ``dpi=`` / ``bbox_inches=`` keep working;
    a ``metadata`` dict may be supplied and is merged with (never clobbered by)
    the channel chunk.
    """
    note = statement(channel, provenance)
    txt = fig.text(0.5, 0.002, footer(channel, provenance),
                   ha="center", va="bottom", fontsize=6.4, color="#b00020",
                   zorder=60, bbox=dict(fc="white", ec="#e0a0a0", lw=0.4,
                                        alpha=0.88))
    fp = _cjk_font()
    if fp is not None:
        txt.set_fontproperties(fp)
    md = dict(kw.pop("metadata", None) or {})
    md[CHUNK] = note
    fig.savefig(path, metadata=md, **kw)
    return note


def read_text(path):
    """PNG tEXt/iTXt chunks of ``path`` as ``{key: value}`` (pure python)."""
    out = {}
    with open(path, "rb") as fh:
        if fh.read(8) != b"\x89PNG\r\n\x1a\n":
            return out                     # not a PNG (e.g. exportfigure JPEG bytes)
        while True:
            hdr = fh.read(8)
            if len(hdr) < 8:
                break
            (n,) = struct.unpack(">I", hdr[:4])
            typ = hdr[4:8]
            data = fh.read(n)
            fh.read(4)                     # CRC
            if typ == b"tEXt":
                k, _, v = data.partition(b"\x00")
                out[k.decode("latin-1")] = v.decode("latin-1", "replace")
            elif typ == b"iTXt":
                key, rest = data.split(b"\x00", 1)
                parts = rest.split(b"\x00", 2)
                if len(parts) == 3:
                    out[key.decode("latin-1")] = parts[2].decode("utf-8", "replace")
            elif typ == b"IEND":
                break
    return out
