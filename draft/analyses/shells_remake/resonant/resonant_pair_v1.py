"""Side-by-side resonant-shell pages (Chris, 11 Oct 22:10): for each scenario, left = equal mass per shell,
right = mass rising outward (x 3^(1/23) per shell, outer = 3x inner, same total). Same widths, same shared log colour scale.
Physics identical to resonant_shells.py (imported functions re-implemented here so that script's outputs aren't rewritten):
N = 24, r_in = 3, R = 15, w_i = c (1 + A psi(x_mid)), A = 0.5, densities / innermost shell of the unmodulated equal-mass set.
Graphs under each wedge (equal height): density (log, same y-range both columns), mass per shell, standing-wave width,
density fit x wave (envelope normalised to the EQUAL-mass fit at r_in in both columns, so the amplitude difference shows).
The plain one-cycle wave graph is dropped: it depends only on the widths, so it is identical in both columns.
Nothing in main.tex or any PDF is touched."""
import os,sys,numpy as np,matplotlib,matplotlib.ticker as mt
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon,Patch
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm,LinearSegmentedColormap
from scipy.optimize import brentq
from scipy.integrate import cumulative_trapezoid
plt.rcParams.update({"font.size":9,"axes.labelsize":8.8,"xtick.labelsize":8,"ytick.labelsize":8})
N=24; R0=3.0; R=15.0; A=0.5; SPAN=45.0; M=1.0; QM=3.0**(1/(N-1))
MI_MO=QM**np.arange(N); MI_MO=N*MI_MO/MI_MO.sum()
CMAP=LinearSegmentedColormap.from_list("shells",["#08154f","#1747a6","#4f86cf","#8a6db0","#b45a8e","#d45a6a","#9b0f1a"])
NORM=LogNorm(10**-1.6,10**0.3)     # = the shared scale of resonant_shells.py (checked below)
SC=[("half_up","Half wavelength going up",r"$\psi=-\cos(\pi x)$",lambda x:-np.cos(np.pi*x)),
    ("half_down","Half wavelength going down",r"$\psi=\cos(\pi x)$",lambda x:np.cos(np.pi*x)),
    ("full_high","Full wavelength, starting high",r"$\psi=\cos(2\pi x)$",lambda x:np.cos(2*np.pi*x)),
    ("full_low","Full wavelength, starting low",r"$\psi=-\cos(2\pi x)$",lambda x:-np.cos(2*np.pi*x))]
SINE=os.environ.get("MODE","cos")=="sine"
if SINE:   # Chris 11 Oct 22:13: particle-in-a-box sine modes, zero at the first and last shell edges
    SC=[("half_pos","n = 1 (+): half wavelength, widest in the middle",r"$\psi=\sin(\pi x)$",lambda x:np.sin(np.pi*x)),
        ("half_neg","n = 1 (−): half wavelength, narrowest in the middle",r"$\psi=-\sin(\pi x)$",lambda x:-np.sin(np.pi*x)),
        ("full_pos","n = 2 (+): full wavelength, widens first",r"$\psi=\sin(2\pi x)$",lambda x:np.sin(2*np.pi*x)),
        ("full_neg","n = 2 (−): full wavelength, narrows first",r"$\psi=-\sin(2\pi x)$",lambda x:-np.sin(2*np.pi*x))]
X=lambda r:(r-R0)/(R-R0)
def march(c,psi):
    e=[R0]
    for i in range(N):
        w=c
        for _ in range(200):
            wn=c*(1+A*psi(X(e[-1]+w/2)))
            if abs(wn-w)<1e-14: break
            w=wn
        e.append(e[-1]+w)
    return np.array(e)
def vol(e): return 4/3*np.pi*(e[1:]**3-e[:-1]**3)
RHO_REF0=M/vol(np.linspace(R0,R,N+1))[0]
def build(psi):
    c=brentq(lambda c: march(c,psi)[-1]-R,0.05,2.0,xtol=1e-14); e=march(c,psi)
    rg=np.linspace(R0,R,20001); ng=cumulative_trapezoid(1/(c*(1+A*psi(X(rg)))),rg,initial=0)
    def fit(r,mo):
        m=MI_MO[0]*QM**(np.interp(r,rg,ng)-0.5) if mo else M
        return m/c/(4*np.pi*r**2*(1+A*psi(X(r))))/RHO_REF0
    return dict(c=c,e=e,w=np.diff(e),fit=fit,rho_eq=M/vol(e)/RHO_REF0,rho_mo=MI_MO/vol(e)/RHO_REF0)
B={k:build(p) for k,_,_,p in SC}
allr=np.concatenate([np.r_[b["rho_eq"],b["rho_mo"]] for b in B.values()]); print("density span",allr.min(),allr.max()); assert NORM.vmin<=allr.min() and allr.max()<=NORM.vmax
# ---------------- geometry (inches) ----------------
W=7.4; LM=0.72; GAPC=0.42; AW=(W-LM-GAPC-0.2)/2; X0=[LM,LM+AW+GAPC]; HW=AW*np.sin(np.radians(SPAN))*1.01
HG=0.62; GP=0.07; NG=4; TOP=0.72; CT=0.26; XL=0.34; LEG=0.62; CBH=0.85
H=TOP+CT+HW+NG*HG+(NG-1)*GP+XL+LEG+CBH
TH=np.linspace(0,np.radians(SPAN),400)
def band(r0,r1): return np.r_[np.c_[r1*np.cos(TH),r1*np.sin(TH)],np.c_[r0*np.cos(TH[::-1]),r0*np.sin(TH[::-1])]]
def waves(e,r):
    i=np.clip(np.searchsorted(e,r,side="right")-1,0,N-1); phi=2*np.pi*(i+(r-e[i])/(e[i+1]-e[i])); return np.where(r>=e[0],np.cos(phi),np.nan)
def page(key,title,lab,psi):
    b=B[key]; e=b["e"]; xe=e/R; w=b["w"]; c=b["c"]
    fig=plt.figure(figsize=(W,H)); Ax=lambda x,y,ww,h: fig.add_axes([x/W,y/H,ww/W,h/H])
    fig.text(0.5,1-0.22/H,(f"Sine mode {title}" if SINE else f"Resonant shells — {title}"),ha="center",va="center",fontsize=12,fontweight="bold")
    fig.text(0.5,1-0.48/H,f"equal mass (left) vs mass rising outward (right);  {N} shells, same widths w ∝ 1 + Aψ(x),  {lab},  A = {A:g}",ha="center",va="center",fontsize=8.8,color="0.25")
    r=np.linspace(0,R,4000); x=r/R; ins=r>=R0; wv=waves(e,r); ref=b["fit"](R0,False)
    for col,(mo,ctitle) in enumerate([(False,"(a) equal mass per shell"),(True,f"(b) mass rising outward (×{QM:.4f}/shell)")]):
        x0=X0[col]; rho=b["rho_mo"] if mo else b["rho_eq"]; mass=MI_MO if mo else np.full(N,M)
        yw=H-TOP-CT-HW; aw=Ax(x0,yw,AW,HW); aw.set_xlim(0,1); aw.set_ylim(0,HW/AW); aw.axis("off")
        aw.set_title(ctitle,loc="left",fontsize=9.2,fontweight="bold",pad=3)
        aw.add_patch(Polygon(np.r_[[[0,0]],np.c_[xe[0]*np.cos(TH),xe[0]*np.sin(TH)]],fc="#dcdcdc",ec="none"))
        for k in range(N): aw.add_patch(Polygon(band(xe[k],xe[k+1]),fc=CMAP(NORM(rho[k])),ec="white",lw=0.45))
        aw.text(0.0,HW/AW*0.99,f"densest {rho.max():.3g}\nemptiest {rho.min():.3g}\nratio {rho.max()/rho.min():.0f}×",fontsize=8.5,va="top",ha="left")
        axs=[]; y=yw
        for _ in range(NG): y-=HG; axs.append(Ax(x0,y,AW,HG)); y-=GP
        a1,a2,a3,a4=axs
        if mo:
            q=b["rho_eq"]; a1.step(np.r_[xe[:-1],xe[-1]],np.r_[q,q[-1]],where="post",color="0.45",lw=0.8,ls="--",alpha=0.6,zorder=2)
        for k in range(N): a1.hlines(rho[k],xe[k],xe[k+1],color=CMAP(NORM(rho[k])),lw=2.4,zorder=3)
        a1.plot(x,np.where(ins,b["fit"](np.maximum(r,1e-9),mo),np.nan),"--",color="0.15",lw=0.8,zorder=4)
        a1.set_yscale("log"); a1.set_ylim(NORM.vmin/1.6,NORM.vmax*1.6); a1.set_yticks([0.03,0.1,0.3,1]); a1.set_yticklabels(["0.03","0.1","0.3","1"]); a1.yaxis.set_minor_formatter(mt.NullFormatter())
        a2.step(np.r_[xe[:-1],xe[-1]],np.r_[mass,mass[-1]],where="post",color="#1e8449",lw=1.4); a2.axhline(1,color="0.6",lw=0.6,ls="--")
        a2.set_ylim(0,2); a2.set_yticks([0,1,2])
        for k in range(N): a3.fill_between([xe[k],xe[k+1]],0,w[k]/c,color=CMAP(NORM(rho[k])),alpha=0.55,lw=0)
        for k in range(N): a3.hlines(w[k]/c,xe[k],xe[k+1],color="0.25",lw=0.9)
        a3.plot(x,np.where(ins,1+A*psi(X(r)),np.nan),color="#6c3483",lw=1.2); a3.axhline(1,color="0.6",lw=0.6); a3.set_ylim(0,1.75); a3.set_yticks([0.5,1,1.5])
        env=np.where(ins,b["fit"](np.maximum(r,1e-9),mo),np.nan)/ref
        a4.plot(x,env,"--",color="#b03050",lw=0.8); a4.plot(x,-env,"--",color="#b03050",lw=0.8); a4.plot(x,env*wv,color="#2c3e50",lw=0.7); a4.axhline(0,color="0.6",lw=0.6)
        a4.set_ylim(-1.3,1.3); a4.set_yticks([-1,0,1])
        if col==0:
            a1.set_ylabel("density\n(÷ ref.)"); a2.set_ylabel("mass per\nshell"); a3.set_ylabel("width ÷ c\n= 1 + Aψ"); a4.set_ylabel("density fit\n× wave")
        for a in axs:
            for bb in xe: a.axvline(bb,color="0.82",lw=0.5,zorder=0)
            a.set_xlim(0,1); a.grid(axis="y",alpha=0.25)
        for a in axs[:-1]: plt.setp(a.get_xticklabels(),visible=False)
        a4.set_xlabel("normalised radius  r / R",labelpad=2)
        ylow=y+GP
    hh=[Line2D([],[],color="0.3",lw=2.4),Line2D([],[],color="0.15",lw=0.8,ls="--"),Line2D([],[],color="0.45",lw=0.8,ls="--",alpha=0.6),
        Line2D([],[],color="#1e8449",lw=1.4),Line2D([],[],color="#6c3483",lw=1.2),Line2D([],[],color="0.25",lw=0.9),
        Line2D([],[],color="#2c3e50",lw=0.7),Line2D([],[],color="#b03050",lw=0.8,ls="--"),Line2D([],[],color="0.82",lw=0.5)]
    ll=["shell density","exact fit  ρ ∝ m(r) / [r²(1 + Aψ)]","equal-mass density (in b)",
        "mass per shell ÷ equal-mass value","standing wave 1 + Aψ(x)","shell widths (steps)",
        "fit × one-cycle-per-shell wave","± fit (÷ equal-mass fit at r$_{in}$)","shell boundaries"]
    fig.legend(hh,ll,loc="upper center",bbox_to_anchor=(0.5,(ylow-XL-0.02)/H),ncol=3,frameon=False,fontsize=8,borderaxespad=0,columnspacing=1.0,handlelength=1.8,labelspacing=0.3)
    cax=fig.add_axes([X0[0]/W,0.45/H,(X0[1]+AW-X0[0])/W,0.14/H])
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=NORM,cmap=CMAP),cax=cax,orientation="horizontal")
    cb.set_label("shell density ÷ innermost shell of the unmodulated equal-mass set  (same log scale as all resonant pages)",fontsize=8.5,labelpad=2); cb.ax.tick_params(labelsize=8)
    fig.savefig(f"resonant_{'sine_' if SINE else ''}pair_{key}.png",dpi=250); plt.close(fig)
for s in SC: page(*s)
print("page in",W,round(H,2))
for k,b in B.items(): print(k,"eq %.3g-%.3g"%(b["rho_eq"].max(),b["rho_eq"].min()),"mo %.3g-%.3g"%(b["rho_mo"].max(),b["rho_mo"].min()))
