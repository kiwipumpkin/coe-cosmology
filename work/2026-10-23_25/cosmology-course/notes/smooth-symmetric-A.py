"""Smooth, axially symmetric A/A_all from our shell geometry (no sine, no boxcar),
at a low and a high observation distance R. Angle psi measured from the anti-CoE
direction (psi=180 points at the CoE), so A(psi)=A(360-psi). Fine 0.25 deg grid."""
import json, os, sys, numpy as np
BUB = "/workspace/cosmology/papers/spreadsheet_model/bubble"
sys.path.insert(0, BUB); cwd = os.getcwd(); os.chdir(BUB)
import ptshell, prof_core as P, bubble_core as BC
os.chdir(cwd)
FINE = np.arange(0, 360, 0.25)
P.TH1 = FINE; BC.TH1 = FINE
fit = json.load(open("migkas-geometry-fit.json"))["wls_fit"]
D, mcc, Rin, Rout = fit["D"], fit["mcc"], fit["Rin"], fit["Rout"]
def sky(R):
    os.chdir(BUB)
    a = dict(DR=D/R, Os=0.049, mcc=float(mcc), rin=Rin/D, rout=Rout/D, Lc=0.005, mc=0.0, Ob=0.0)
    rr = P.run("ptshell", a, 0.0, 0.0); os.chdir(cwd)
    af = rr["af"]; return af/af.mean(), rr["vrmin"]
out = {}
for R in (3.0, 14.30048287260589):
    A, vr = sky(R); out[R] = A
    sym = np.max(np.abs(A - np.roll(A[::-1], 1)))
    i = np.argmax(A); print(f"R={R:.2f}  D={D:.2f} Rin={Rin:.2f} Rout={Rout:.2f}  vrmin={vr:.3g}  sym_err={sym:.2e}")
    # extrema
    n=len(A)
    for k in range(n):
        a0,a1,a2=A[k-1],A[k],A[(k+1)%n]
        if a1>a0 and a1>=a2: print("  max", FINE[k], round(a1,4))
        if a1<a0 and a1<=a2: print("  min", FINE[k], round(a1,4))
np.savez("smooth-symmetric-A.npz", psi=FINE, **{f"R{R:.1f}": v for R, v in out.items()})
import csv
with open("smooth-symmetric-A.csv","w",newline="") as f:
    w=csv.writer(f); keys=list(out)
    w.writerow(["psi_deg"]+[f"A_over_Aall_R{R:.1f}Gpc" for R in keys])
    for j,p in enumerate(FINE): w.writerow([f"{p:.2f}"]+[f"{out[R][j]:.6f}" for R in keys])
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":15,"mathtext.fontset":"cm"})
fig,ax=plt.subplots(figsize=(11,6))
for R,c in zip(out,("C0","C3")):
    ax.plot(np.append(FINE,360), np.append(out[R],out[R][0]), c, lw=2.4,
            label=rf"$R$ = {R:.1f} Gpc" + (" (low)" if R<10 else " (our fit)"))
ax.axvline(180,color="0.5",ls=":",lw=1.2)
ax.set_xlim(0,360); ax.set_xticks(np.arange(0,361,45))
ax.set_xlabel(r"angle from anti-CoE direction $\psi$ (deg); $\psi=180^\circ$ points at the CoE")
ax.set_ylabel(r"$A/A_{\rm all}$"); ax.grid(alpha=0.3)
ax.legend(loc="lower center",bbox_to_anchor=(0.5,1.0),ncol=2,frameon=False,fontsize=15)
ax.set_title(rf"Shell geometry, no sine: $D$={D:.1f} Gpc, shell {Rin:.1f}–{Rout:.1f} Gpc",pad=48,fontsize=15)
fig.tight_layout(); fig.savefig("smooth-symmetric-A.png",dpi=200)
print("ok")
