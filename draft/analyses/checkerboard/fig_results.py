import json,numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import sys; sys.path.insert(0,"."); 
d=np.load("final.npz"); F=float(d["FLOOR"]); R=d["R"]; res=json.load(open("final.json")); ens=json.load(open("runs.json"))
AXIS=15.16; l=np.arange(360)
def modes(w):
    th=np.radians(l-AXIS); m=np.arange(1,181)
    a=2*np.mean(w[None]*np.cos(m[:,None]*th),1); b=2*np.mean(w[None]*np.sin(m[:,None]*th),1)
    return m,np.sqrt((a**2+b**2)/2)/F,np.sqrt(a**2/2)/F
def lowp(w,M=12):
    th=np.radians(l-AXIS); m=np.arange(1,M+1)
    a=2*np.mean(w[None]*np.cos(m[:,None]*th),1); b=2*np.mean(w[None]*np.sin(m[:,None]*th),1)
    return (a[:,None]*np.cos(m[:,None]*th)+b[:,None]*np.sin(m[:,None]*th)).sum(0)
COL={"chess":"#d62728","einasto":"#7b3294"}; NAME={"chess":"chessboard (85.7 Mpc cells)","einasto":"Einasto cell (void centre, clusters on walls, 171 Mpc)"}
def band(ax): ax.axhspan(-1,1,color="#2ca02c",alpha=0.28,ec="#2ca02c",lw=0.8,label="±1× (CMB level, 10⁻⁵)")
def leg(ax,y=-0.2,n=2): ax.legend(loc="upper center",bbox_to_anchor=(0.5,y),ncol=n,fontsize=8.5,frameon=False)
def mark(ax): 
    for x in (15.16,195.16): ax.axvline(x,color="k",ls=":",lw=0.8)
# ---------------- figure 4: checkerboard alone
fig,axs=plt.subplots(3,1,figsize=(10,16),gridspec_kw=dict(hspace=0.68))
ax=axs[0]; band(ax)
for k in ["chess","einasto"]: ax.plot(l,d[f"{k}_w"]/F,color=COL[k],lw=0.9,label=f"{NAME[k]}: {res[k]['own']:.1f}× rms")
mark(ax); ax.set_xlim(0,360); ax.set_xlabel("ring position l (deg); dotted = CoE axis (15°/195°)"); ax.set_ylabel("ring residual (× 10⁻⁵)")
ax.set_title("(a) CMB ring residual from the periodic background alone (contrast A = 0.10, Einasto/supergalactic orientation)",fontsize=10); leg(ax,-0.2,1)
ax=axs[1]
for k,off in [("chess",-0.2),("einasto",0.2)]:
    m,amp,_=modes(d[f"{k}_w"]); ax.bar(m+off,amp,width=0.4,color=COL[k],label=f"{NAME[k]}: m1–12 {res[k]['own_m1_12']:.1f}×, m13+ {res[k]['own_13']:.1f}×, symmetric share {res[k]['own_sym']*100:.0f}%")
ax.axhline(1,color="#2ca02c",lw=2,label="1× (CMB level)"); ax.axvspan(0.5,12.5,color="0.85",alpha=0.5,label="m = 1–12 (the modes we solve for)")
ax.set_xlim(0,181); ax.set_yscale("log"); ax.set_ylim(1e-2,10); ax.set_xlabel("ring mode m (1° sampling, so m ≤ 180)"); ax.set_ylabel("amplitude per mode (×)")
ax.set_title("(b) Spectrum by mode: the cells themselves map to m ≈ 524, so their power aliases into a flat, noise-like floor at all m;\nthe chessboard's m=1 bar is mostly the pull of the few cells right around us (a dipole)",fontsize=10); leg(ax,-0.2,1)
ax=axs[2]; Aa=np.array([0.01,0.03,0.10,0.19,0.3])
for k in ["chess","einasto"]:
    ax.plot(Aa,res[k]["own_13"]*Aa/0.10,"-o",color=COL[k],label=f"{NAME[k]}: m13+ noise")
    ax.plot(Aa,res[k]["own_m1_12"]*Aa/0.10,"--s",color=COL[k],mfc="w",label=f"{NAME[k]}: m1–12")
ax.axhline(1,color="#2ca02c",lw=2,label="1× (CMB level)"); ax.axvline(0.10,color="k",ls=":",lw=1,label="A = 0.10 used (≈ periodic share of σ = 0.19 at 37 h⁻¹ Mpc)")
ax.axvline(0.19,color="0.5",ls=":",lw=1,label="A = 0.19 = all fluctuations at the cell scale")
ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel("cell density contrast A (rms)"); ax.set_ylabel("ring residual (×)")
ax.set_title("(c) Scaling: the pattern is linear in A (and in the matter density used; baryons only would be ~6× lower)",fontsize=10); leg(ax,-0.2,2)
fig.savefig("fig4_checkerboard_ring_alone.png",dpi=130,bbox_inches="tight"); plt.close(fig)
# ---------------- figure 5: add bodies, show the noise left
fig,axs=plt.subplots(4,1,figsize=(10,23),gridspec_kw=dict(hspace=0.78,height_ratios=[1,1,1,1.25]))
ax=axs[0]; band(ax); ax.plot(l,R/F,"k-",lw=1.5,label=f"SR + central mass residual alone: {np.sqrt(np.mean(R**2))/F:.0f}×")
for k in ["chess","einasto"]: ax.plot(l,(R+d[f'{k}_w'])/F,color=COL[k],lw=0.8,label=f"+ {NAME[k]}: {res[k]['before']:.0f}×")
mark(ax); ax.set_xlim(0,360); ax.set_ylabel("ring residual (× 10⁻⁵)"); ax.set_xlabel("ring position l (deg)")
ax.set_title("(a) Before: checkerboard on top of the 237× residual (the checkerboard is the thin jitter on the big wave)",fontsize=10); leg(ax,-0.2,1)
ax=axs[1]; band(ax)
for k in ["chess","einasto"]:
    t=d[f"{k}_tot"]; ax.plot(l,t/F,color=COL[k],lw=0.8,label=f"after bodies, {NAME[k]}: total {res[k]['after']:.1f}×")
    ax.plot(l,lowp(t)/F,color=COL[k],lw=2.5,alpha=0.9,ls="--",label=f"   its m = 1–12 part: {res[k]['after_m1_12']:.3f}×")
mark(ax); ax.set_xlim(0,360); ax.set_ylabel("ring residual (× 10⁻⁵)"); ax.set_xlabel("ring position l (deg)")
ax.set_title(f"(b) After: {res['chess']['n_bodies']} bodies (≥ {res['chess']['min_ring_dist_Gpc']:.1f} Gpc from every ring point) cancel m = 1–12; the checkerboard noise stays",fontsize=10); leg(ax,-0.2,1)
ax=axs[2]
for k,off in [("chess",-0.2),("einasto",0.2)]:
    m,amp,_=modes(d[f"{k}_tot"]); ax.bar(m+off,amp,width=0.4,color=COL[k],label=f"{NAME[k]}: m1–12 {res[k]['after_m1_12']:.3f}×, m13+ {res[k]['after_13']:.1f}×, symmetric share {res[k]['after_sym']*100:.0f}%")
ax.axhline(1,color="#2ca02c",lw=2,label="1× (CMB level)"); ax.axvspan(0.5,12.5,color="0.85",alpha=0.5,label="m = 1–12: cancelled")
ax.set_xlim(0,181); ax.set_yscale("log"); ax.set_ylim(1e-4,10); ax.set_xlabel("ring mode m"); ax.set_ylabel("amplitude per mode (×)")
ax.set_title("(c) Spectrum after bodies: m ≤ 12 pushed far below 1×; m ≥ 13 is the checkerboard's flat, noise-like floor",fontsize=10); leg(ax,-0.2,1)
ax=axs[3]; X=d["chess_X"]; M=d["chess_M"]; s=np.abs(M)>1e47
th=np.linspace(0,2*np.pi,400); ax.fill_between(np.cos(th)*15.8,np.sin(th)*15.8,np.sin(th)*15.8*0,color="none")
from matplotlib.patches import Circle, Wedge
ax.add_patch(Wedge((0,0),14.3+1.5,0,360,width=3.0,color="#2ca02c",alpha=0.3,label="kept clear: within 1.5 Gpc of the ring"))
ax.add_patch(Circle((0,0),14.3,fill=False,color="#1a9e3f",lw=2,label="CMB ring (14.3 Gpc)"))
ax.add_patch(Circle(tuple(d.get("C",np.array([-56.17,-15.22]))[:2]) if False else (-56.17,-15.22),43.9,fill=False,color="m",ls="--",lw=1,label="CoE core edge (43.9 Gpc from CoE)"))
pos=s&(M>0); neg=s&(M<0)
ax.scatter(X[pos,0],X[pos,1],s=8+40*np.abs(M[pos])/2e52,c="#d62728",alpha=0.7,label=f"added masses ({pos.sum()})")
ax.scatter(X[neg,0],X[neg,1],s=8+40*np.abs(M[neg])/2e52,c="#3b6fd8",alpha=0.7,label=f"added voids / negative masses ({neg.sum()})")
ax.plot(0,0,"k*",ms=12,label="us"); ax.annotate("to CoE",xy=(-29,-7.9),xytext=(-18,-4.9),arrowprops=dict(arrowstyle="-|>"))
ax.set_xlim(-31,31); ax.set_ylim(-31,31); ax.set_aspect("equal"); ax.set_xlabel("Gpc (ring plane)"); ax.set_ylabel("Gpc")
ax.set_title("(d) Side view of the added bodies, chessboard case (projected onto the ring plane: dots over the green band sit\nabove or below it; in 3D every body is ≥ 4.7 Gpc from the ring). Dot size ∝ |mass|, max 2×10⁵² kg",fontsize=10); leg(ax,-0.1,2)
fig.savefig("fig5_bodies_cancel_noise_left.png",dpi=130,bbox_inches="tight")
