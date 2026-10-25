"""Log-scale versions of the extrema overview. Same tracked extrema as
extrema-vs-inputs.py (imported unchanged); only the axes change."""
import importlib.util, sys
import numpy as np, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
NOTES = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("ex", NOTES / "extrema-vs-inputs.py")
ex = importlib.util.module_from_spec(spec); sys.modules["ex"] = ex; spec.loader.exec_module(ex)
plt.rcParams.update({"mathtext.fontset": "dejavusans", "font.size": 14})

Ts = [ex.track(key, lo, hi, base, over) for key, label, lo, hi, base, tag, over in ex.SWEEPS]

def make(logx, fname):
    fig, axs = plt.subplots(2, 4, figsize=(16, 8.6), sharex="col",
                            gridspec_kw=dict(hspace=0.08, wspace=0.30))
    for c, (T, (key, label, lo, hi, base, tag, over)) in enumerate(zip(Ts, ex.SWEEPS)):
        al, aa = axs[0, c], axs[1, c]
        for b in T["branches"]:
            ls = "-" if b["kind"] == "max" else "--"
            al.plot(T["vals"], b["l"] % 360, ls=ls, color=b["col"], lw=2.6)
            aa.plot(T["vals"], np.abs(b["A"] - 1.0), ls=ls, color=b["col"], lw=2.6)
        for m in T["merges"]:
            al.plot(m["value"], m["l_deg"], "o", mfc="none", mec="k", mew=2, ms=11, zorder=6)
            aa.plot(m["value"], abs(m["A"] - 1), "o", mfc="none", mec="k", mew=2, ms=11, zorder=6)
        for a in (al, aa):
            a.axvline(base, color="k", ls=":", lw=1.8)
            a.set_yscale("log"); a.grid(alpha=0.35, which="both")
            if logx:
                from matplotlib.ticker import FixedLocator, NullLocator, FormatStrFormatter
                a.set_xscale("log")
                tk = {"x":[3,4,6,8,10],"alpha":[0.2,0.5,1,2],"kappa":[0.15,0.2,0.3,0.5,0.7],"eps":[0.06,0.1,0.2,0.4,0.7]}[key]
                a.xaxis.set_major_locator(FixedLocator(tk)); a.xaxis.set_minor_locator(NullLocator())
                a.xaxis.set_major_formatter(FormatStrFormatter("%g"))
            a.set_xlim(T["vals"][0], T["vals"][-1]); a.tick_params(labelsize=13)
        al.set_ylim(2, 400); aa.set_ylim(3e-3, 0.6)
        al.set_title(f"vary {label}" + (r" ($\alpha$=1.3 held)" if key == "eps" else ""), fontsize=15)
        aa.set_xlabel(label + (" (log)" if logx else ""), fontsize=15)
    axs[0, 0].set_ylabel(r"$\ell$ of extremum [deg] (log)", fontsize=15)
    axs[1, 0].set_ylabel(r"$|A/A_{\mathrm{all}}-1|$ at extremum (log)", fontsize=15)
    hs = [Line2D([], [], color=ex.COLS["max"][0], lw=2.6, label="max near 10°"),
          Line2D([], [], color=ex.COLS["max"][1], lw=2.6, label="max near 170°"),
          Line2D([], [], color=ex.COLS["min"][0], lw=2.6, ls="--", label="min near 95°"),
          Line2D([], [], color=ex.COLS["min"][1], lw=2.6, ls="--", label="min near 266°"),
          Line2D([], [], color="k", lw=1.8, ls=":", label="baseline"),
          Line2D([], [], marker="o", mfc="none", mec="k", mew=2, ms=11, ls="none", label="max–min merge")]
    fig.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, 0.995), ncol=6,
               fontsize=14, frameon=False, columnspacing=1.2, handlelength=2.2)
    fig.subplots_adjust(top=0.86, bottom=0.10, left=0.07, right=0.99)
    fig.savefig(NOTES / fname, dpi=200); plt.close(fig)

make(False, "extrema-vs-inputs-logy.png")
make(True, "extrema-vs-inputs-loglog.png")
print("ok")
