"""Standalone figure: radial 150 Mpc wave and 120 Mpc quasi-lattice backgrounds alone (no bodies added), side by side."""
import numpy as np,json,matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle,Rectangle
from matplotlib.lines import Line2D
from common import *
import sys; sys.argv=["x","none_s0"]
from fit_on_web import spec
plt.rcParams.update({"font.size":7.5})
def field(kind,seed):     # same construction / random draws as wave_bg.py
    rng=np.random.default_rng(seed); A=0.15; LAM=0.15 if kind=="wave" else 0.12; k=2*np.pi/LAM; JIT=0.25
    NW=40; kv=rng.normal(size=(NW,3)); kv/=np.linalg.norm(kv,axis=1,keepdims=True); kv*=(2*np.pi/rng.uniform(0.5,2.0,NW))[:,None]; ph=rng.uniform(0,2*np.pi,NW)
    jit=lambda P,c=0: np.sqrt(2/NW)*np.cos(P@np.roll(kv,c,0).T+np.roll(ph,7*c)).sum(1)
    Q,_=np.linalg.qr(rng.normal(size=(3,3)))
    def delta(P):
        if kind=="wave": r=np.linalg.norm(P-C_G,axis=1); return A*np.cos(k*r+2*np.pi*JIT*jit(P))
        Pp=P@Q+LAM*JIT*np.c_[jit(P,0),jit(P,1),jit(P,2)]; return A/3*(np.cos(k*Pp[:,0])+np.cos(k*Pp[:,1])+np.cos(k*Pp[:,2]))
    return delta
def img(delta,c,h,n):
    x=np.linspace(-h,h,n); X,Y=np.meshgrid(x,x,indexing="ij"); P=np.c_[X.ravel()+c[0],Y.ravel()+c[1],np.zeros(X.size)]
    out=np.zeros(len(P))
    for i in range(0,len(P),200000): out[i:i+200000]=delta(P[i:i+200000])
    return out.reshape(X.shape).T
COL={"wave":"#d35400","lattice":"#16a085"}; NM={"wave":"Radial wave: concentric shells about the CoE, 150 Mpc wavelength","lattice":"Jittered quasi-lattice, 120 Mpc spacing"}
c195=np.array([np.cos(np.radians(195)),np.sin(np.radians(195)),0])*RMIN
fig=plt.figure(figsize=(10.5,15)); G=fig.add_gridspec(5,2,height_ratios=[1.25,1.1,1.1,0.7,0.9],hspace=0.38,wspace=0.18,left=0.07,right=0.98,top=0.95,bottom=0.07)
th=np.linspace(0,2*np.pi,800)
for j,kind in enumerate(("wave","lattice")):
    d=field(kind,1); col=COL[kind]; W=[np.load(f"{OUT}/final_{kind}_s{s}.npz")["w"] for s in range(1,6)]; S=[spec(w) for w in W]
    # overview
    ax=fig.add_subplot(G[0,j]); ax.add_patch(Circle(C_G[:2],A_IN_G,fc="#fdf2e9",ec="#e6b0aa",ls="--",lw=0.6))
    if kind=="wave":
        rr=np.linalg.norm(C_G[:2]); 
        for q in np.arange(0.15*5,90,0.15*40): ax.plot(C_G[0]+q*np.cos(th),C_G[1]+q*np.sin(th),color=col,lw=0.35,alpha=0.7)
        ax.text(-104,43,"every 40th crest drawn (actual crests every 150 Mpc)\nshells fill the core, the shell and our sphere",fontsize=6,color=col,va="top",bbox=dict(fc="white",ec="none",alpha=0.85,pad=1))
    else:
        ax.add_patch(Rectangle((-106,-64),154,112,fc=col,alpha=0.08,ec="none"))
        ax.text(-104,43,"lattice fills the whole region (≈10$^{9}$ nodes,\n120 Mpc apart — too fine to draw at this scale)",fontsize=6,color=col,va="top",bbox=dict(fc="white",ec="none",alpha=0.85,pad=1))
    ax.add_patch(Circle((0,0),RMIN,fill=False,ec="#1f4e79",lw=1.3)); ax.plot(0,0,"k+")
    ax.add_patch(Rectangle(c195[:2]-6,12,12,fill=False,ec="k",lw=0.8)); ax.text(c195[0]-6,c195[1]+6.5,"zoom A",fontsize=6)
    ax.set_xlim(-106,48); ax.set_ylim(-64,48); ax.set_aspect("equal"); ax.set_title(NM[kind],color=col,fontweight="bold",fontsize=8); ax.set_xlabel("Gpc (ring plane, seen face-on)")
    ax.text(C_G[0],C_G[1]-A_IN_G+3,"CoE core",fontsize=6,ha="center",color="#a93226")
    for s_ in ax.spines.values(): s_.set_edgecolor(col); s_.set_linewidth(1.8)
    # zoom A (12 Gpc) and zoom B (1.5 Gpc)
    for row,h,n,lab in ((1,6.0,1500,"zoom A: 12 Gpc box on the ring at ℓ=195°"),(2,0.75,600,"zoom B: 1.5 Gpc box at the same point")):
        ax=fig.add_subplot(G[row,j]); I=img(d,c195,h,n)
        ax.imshow(I,origin="lower",extent=[c195[0]-h,c195[0]+h,c195[1]-h,c195[1]+h],cmap="RdBu_r",vmin=-0.15,vmax=0.15,interpolation="bilinear")
        ax.plot(RMIN*np.cos(th),RMIN*np.sin(th),color="k",lw=1.4); ax.set_xlim(c195[0]-h,c195[0]+h); ax.set_ylim(c195[1]-h,c195[1]+h); ax.set_aspect("equal")
        ax.set_title(lab+" (black = CMB ring)",fontsize=7); ax.tick_params(labelsize=6)
        if row==2: ax.add_patch(Rectangle((c195[0]+h-0.42,c195[1]-h+0.06),0.15 if kind=="wave" else 0.12,0.025,fc="k")); ax.text(c195[0]+h-0.42,c195[1]-h+0.1,"150 Mpc" if kind=="wave" else "120 Mpc",fontsize=6)
    # ring pattern of the background alone
    ax=fig.add_subplot(G[3,j]); ax.axhspan(-1,1,color="green",alpha=0.3,lw=0); ax.axhline(0,color="0.6",lw=0.5)
    ax.plot(np.arange(360),W[0]/FLOOR,color=col,lw=0.8); lim=1.15*max(np.abs(W[0]/FLOOR).max(),1); ax.set_ylim(-lim,lim); ax.set_xlim(0,360); ax.set_xticks([0,90,180,270,360])
    ax.set_title(f"its ring pattern alone (seed 1): {S[0]['total']:.1f}× rms; symmetric share of m=1–12: {S[0]['sym']*100:.0f}%",fontsize=7,color=col); ax.set_xlabel("ℓ (°) around the ring"); ax.set_ylabel("× $10^{-5}$")
    # spectrum, sym/anti stacked, median over 5 seeds
    ax=fig.add_subplot(G[4,j]); amp=np.array([s_["amp"] for s_ in S]); sm=np.array([s_["sym_m"] for s_ in S]); m=np.arange(1,13)
    a_med=np.median(amp,0); sym_med=np.median(sm,0)
    ax.bar(m,a_med*sym_med,color=col,width=0.7,label="mirror-symmetric part"); ax.bar(m,a_med*(1-sym_med),bottom=a_med*sym_med,color="0.75",width=0.7,label="antisymmetric part")
    ax.errorbar(m,a_med,yerr=[a_med-amp.min(0),amp.max(0)-a_med],fmt="none",ecolor="k",lw=0.6,capsize=1.5)
    h13=np.array([s_["rms_13up"] for s_ in S]); ax.bar([14],[np.median(h13)],color="0.35",width=0.9,label="m ≥ 13 (rms)"); ax.errorbar([14],[np.median(h13)],yerr=[[np.median(h13)-h13.min()],[h13.max()-np.median(h13)]],fmt="none",ecolor="k",lw=0.6)
    ax.axhline(1,color="green",lw=1.5,label="1× (CMB level)"); ax.set_xticks(list(m)+[14]); ax.set_xticklabels([str(x) for x in m]+["13+"]); ax.set_ylabel("amplitude (× floor)")
    ax.set_title("by-mode spectrum, median of 5 seeds (bars = range)",fontsize=7); ax.set_xlabel("mode m around the ring")
    ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.28),ncol=4,frameon=False,fontsize=6.5)
hd=[Line2D([],[],color="#1f4e79",lw=1.3,label="CMB ring / our sphere (14.3 Gpc)"),Rectangle((0,0),1,1,fc="#fdf2e9",ec="#e6b0aa",ls="--",label="CoE core"),
    Rectangle((0,0),1,1,fc="#b2182b",label="over-dense (δ up to +0.15)"),Rectangle((0,0),1,1,fc="#2166ac",label="under-dense (δ down to −0.15)"),Rectangle((0,0),1,1,fc="green",alpha=0.3,label="±1× (CMB level)")]
fig.legend(handles=hd,loc="lower center",ncol=5,frameon=False,fontsize=7,bbox_to_anchor=(0.5,0.0))
fig.savefig(f"{OUT}/fig_wave_lattice_alone.png",dpi=170,bbox_inches="tight"); fig.savefig(f"{OUT}/fig_wave_lattice_alone.pdf",bbox_inches="tight")
print({k:[ (round(s['total'],2),round(s['sym'],2)) for s in [spec(np.load(f'{OUT}/final_{k}_s{i}.npz')['w']) for i in range(1,6)]] for k in ('wave','lattice')})
