"""Resonant (particle-in-a-box) shell widths, 24 equal-mass shells (Chris, 11 Oct 2026).
Geometry: N = 24 shells between r_in = 3 and R = 15 (same inner radius ratio r_in/R = 0.2 as the constant-spacing figure).
Standing-wave modulation over the shell region, x = (r - r_in)/(R - r_in) in [0,1]:
    half_up    psi = -cos(pi x)     half_down  psi = cos(pi x)
    full_high  psi =  cos(2 pi x)   full_low   psi = -cos(2 pi x)
Shell widths  w_i = c * (1 + A psi(x_mid,i)),  A = 0.5 (widest/narrowest local width = 3),
    x_mid,i evaluated at the shell's own midpoint (implicit, solved by fixed-point per shell),
    c chosen (root-find) so the 24 shells exactly fill r_in..R.
Equal mass per shell: rho_i = M / [(4/3) pi (r_{i+1}^3 - r_i^3)], shown divided by the innermost-shell density of the UNMODULATED
    24-shell reference (A = 0, constant width 0.5), so one colour scale compares all four pages.
Smooth fit: continuum limit of equal mass per shell, rho(r) = (M/c) / [4 pi r^2 (1 + A psi(x))] (no free parameters once c is known).
Wave: one cycle per shell, crests on every boundary: cos(2 pi * shells crossed).
Bottom: rho_fit(r)/rho_fit(r_in) * wave, with +/- envelope.
Nothing in main.tex or any PDF is touched."""
import json,numpy as np,matplotlib,matplotlib.ticker as mt
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.lines import Line2D
from matplotlib.colors import LogNorm,LinearSegmentedColormap
from scipy.optimize import brentq
plt.rcParams.update({"font.size":9.5,"axes.labelsize":9.5,"xtick.labelsize":8.5,"ytick.labelsize":8.5})
N=24; R0=3.0; R=15.0; A=0.5; SPAN=45.0; M=1.0
CMAP=LinearSegmentedColormap.from_list("shells",["#08154f","#1747a6","#4f86cf","#8a6db0","#b45a8e","#d45a6a","#9b0f1a"])
SC=[("half_up","Half wavelength going up",r"$\psi=-\cos(\pi x)$",lambda x:-np.cos(np.pi*x)),
    ("half_down","Half wavelength going down",r"$\psi=\cos(\pi x)$",lambda x:np.cos(np.pi*x)),
    ("full_high","Full wavelength, starting high",r"$\psi=\cos(2\pi x)$",lambda x:np.cos(2*np.pi*x)),
    ("full_low","Full wavelength, starting low",r"$\psi=-\cos(2\pi x)$",lambda x:-np.cos(2*np.pi*x))]
COLS=dict(half_up="#b03050",half_down="#1747a6",full_high="#1e8449",full_low="#d68910")
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
def solve(psi):
    c=brentq(lambda c: march(c,psi)[-1]-R,0.05,2.0,xtol=1e-14); e=march(c,psi); return c,e
def vol(e): return 4/3*np.pi*(e[1:]**3-e[:-1]**3)
eref=np.linspace(R0,R,N+1); RHO_REF0=M/vol(eref)[0]
RES={}
for key,title,lab,psi in SC:
    c,e=solve(psi); rho=M/vol(e)/RHO_REF0; w=np.diff(e); xm=X(0.5*(e[1:]+e[:-1]))
    assert np.allclose(w,c*(1+A*psi(xm)),atol=1e-10) and abs(e[-1]-R)<1e-10
    fit=lambda r,c=c,psi=psi:(M/c)/(4*np.pi*r**2*(1+A*psi(X(r))))/RHO_REF0
    rm=0.5*(e[1:]+e[:-1]); err=np.abs(fit(rm)/rho-1).max()
    RES[key]=dict(c=c,e=e,rho=rho,w=w,fit=fit,psi=psi,title=title,lab=lab,fit_max_rel_err=float(err))
allr=np.concatenate([v["rho"] for v in RES.values()])
NORM=LogNorm(10**(np.floor(np.log10(allr.min())*10)/10),10**(np.ceil(np.log10(allr.max())*10)/10))
print("norm",NORM.vmin,NORM.vmax)
# ---------------- page geometry (inches) ----------------
W=7.4; LM=0.95; AW=6.0; HW=AW*np.sin(np.radians(SPAN))*1.01; HG=0.66; GP=0.08
TOP=0.62; XL=0.36; LEGH=0.62; CB=0.80
H=TOP+HW+4*HG+3*GP+XL+LEGH+CB
TH=np.linspace(0,np.radians(SPAN),400)
def band(r0,r1): return np.r_[np.c_[r1*np.cos(TH),r1*np.sin(TH)],np.c_[r0*np.cos(TH[::-1]),r0*np.sin(TH[::-1])]]
def waves(e,r):
    i=np.clip(np.searchsorted(e,r,side="right")-1,0,N-1); phi=2*np.pi*(i+(r-e[i])/(e[i+1]-e[i])); return np.where(r>=e[0],np.cos(phi),np.nan)
def page(key,out):
    d=RES[key]; e=d["e"]; xe=e/R; rho=d["rho"]; w=d["w"]; c=d["c"]
    fig=plt.figure(figsize=(W,H)); Ax=lambda y,h: fig.add_axes([LM/W,y/H,AW/W,h/H])
    fig.text(0.5,1-0.22/H,f"Resonant shell widths — {d['title']}",ha="center",va="center",fontsize=12.5,fontweight="bold")
    fig.text(0.5,1-0.45/H,f"{N} equal-mass shells; local width ∝ 1 + A ψ(x),  {d['lab']},  A = {A:g},  x = (r − r$_{{in}}$)/(R − r$_{{in}}$)",ha="center",va="center",fontsize=9.5,color="0.25")
    yw=H-TOP-HW; aw=Ax(yw,HW); aw.set_xlim(0,1); aw.set_ylim(0,HW/AW); aw.axis("off")
    aw.add_patch(Polygon(np.r_[[[0,0]],np.c_[xe[0]*np.cos(TH),xe[0]*np.sin(TH)]],fc="#dcdcdc",ec="none"))
    for k in range(N): aw.add_patch(Polygon(band(xe[k],xe[k+1]),fc=CMAP(NORM(rho[k])),ec="white",lw=0.6))
    aw.text(0.0,HW/AW*0.985,f"{N} shells, equal mass each\ninner radius {R0:g}, outer radius R = {R:g}\nwidths {w.min():.3f} – {w.max():.3f}  (mean {w.mean():.3f})\ndensity {rho.max()/rho.min():.0f}× from densest to emptiest shell",fontsize=9.5,va="top",ha="left")
    axs=[]; y=yw
    for _ in range(4): y-=HG; axs.append(Ax(y,HG)); y-=GP
    a1,a2,a3,a4=axs; r=np.linspace(0,R,4000); x=r/R; ins=r>=R0
    # density
    a1.step(np.r_[eref[:-1],eref[-1]]/R,np.r_[M/vol(eref)/RHO_REF0,(M/vol(eref)/RHO_REF0)[-1]],where="post",color="0.55",lw=0.8,ls=":",zorder=2)
    for k in range(N): a1.hlines(rho[k],xe[k],xe[k+1],color=CMAP(NORM(rho[k])),lw=2.8,zorder=3)
    rf=np.where(ins,d["fit"](np.maximum(r,1e-9)),np.nan); a1.plot(x,rf,"--",color="0.15",lw=0.9,zorder=4)
    a1.set_yscale("log"); a1.set_ylim(NORM.vmin/1.6,NORM.vmax*1.6); a1.set_ylabel("density\n(÷ ref. inner)",fontsize=9)
    a1.yaxis.set_major_locator(mt.LogLocator(numticks=4)); a1.yaxis.set_minor_formatter(mt.NullFormatter())
    # width modulation
    rel=w/c
    for k in range(N): a2.fill_between([xe[k],xe[k+1]],0,rel[k],color=CMAP(NORM(rho[k])),alpha=0.55,lw=0,step=None)
    for k in range(N): a2.hlines(rel[k],xe[k],xe[k+1],color="0.25",lw=1.0)
    a2.plot(x,np.where(ins,1+A*d["psi"](X(r)),np.nan),color="#6c3483",lw=1.3)
    a2.axhline(1,color="0.6",lw=0.6); a2.set_ylim(0,1.75); a2.set_yticks([0.5,1,1.5]); a2.set_ylabel("width ÷ c\n= 1 + Aψ",fontsize=9)
    t=a2.twinx(); t.set_ylim((0-1)/A,(1.75-1)/A); t.set_yticks([-1,0,1]); t.set_ylabel("ψ",fontsize=9,rotation=0,labelpad=6); t.tick_params(labelsize=8.5)
    # wave
    wv=waves(e,r); a3.plot(x,wv,color="#2c3e50",lw=0.8); a3.set_ylim(-1.3,1.3); a3.set_yticks([-1,0,1]); a3.set_ylabel("wave",fontsize=9)
    env=rf/d["fit"](R0); a4.plot(x,env,"--",color="#b03050",lw=0.9); a4.plot(x,-env,"--",color="#b03050",lw=0.9); a4.plot(x,env*wv,color="#2c3e50",lw=0.8)
    a4.set_ylim(-1.3,1.3); a4.set_yticks([-1,0,1]); a4.set_ylabel("density fit\n× wave",fontsize=9)
    for a in (a3,a4): a.axhline(0,color="0.6",lw=0.6)
    for a in axs:
        for b in xe: a.axvline(b,color="0.8",lw=0.6,zorder=0)
        a.set_xlim(0,1); a.grid(axis="y",alpha=0.25)
    for a in axs[:-1]: plt.setp(a.get_xticklabels(),visible=False)
    a4.set_xlabel("normalised radius  r / R   (grey lines = shell boundaries)",labelpad=2)
    hh=[Line2D([],[],color="0.3",lw=2.8),Line2D([],[],color="0.15",lw=0.9,ls="--"),Line2D([],[],color="0.55",lw=0.8,ls=":"),
        Line2D([],[],color="#6c3483",lw=1.3),Line2D([],[],color="0.25",lw=1.0),Line2D([],[],color="#2c3e50",lw=0.8),
        Line2D([],[],color="#b03050",lw=0.9,ls="--"),Line2D([],[],color="0.8",lw=0.6)]
    ll=["shell density","smooth fit  ρ ∝ 1 / [r² (1 + Aψ)]","unmodulated 24 shells (A = 0)",
        "standing wave 1 + Aψ(x)","shell widths (steps)","wave (one cycle per shell); fit × wave","± density fit","shell boundaries"]
    fig.legend(hh,ll,loc="upper center",bbox_to_anchor=(0.5,(y+GP-XL-0.04)/H),ncol=3,frameon=False,fontsize=8.8,borderaxespad=0,columnspacing=1.4,handlelength=2.0,labelspacing=0.3)
    cax=fig.add_axes([LM/W,0.45/H,AW/W,0.15/H])
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=NORM,cmap=CMAP),cax=cax,orientation="horizontal")
    cb.set_label("shell density ÷ innermost shell of the unmodulated 24-shell set  (same log scale on all four pages)",fontsize=9,labelpad=2); cb.ax.tick_params(labelsize=8.5)
    fig.savefig(out,dpi=250); plt.close(fig)
for key,*_ in SC: page(key,f"resonant_{key}.png")
# ---------------- comparison ----------------
fig,(b0,b1)=plt.subplots(2,1,figsize=(7.4,8.6),sharex=True,gridspec_kw=dict(hspace=0.08,left=0.12,right=0.97,top=0.90,bottom=0.22))
fig.suptitle(f"Resonant shell widths: the four standing-wave scenarios\n{N} equal-mass shells, width ∝ 1 + Aψ(x), A = {A:g}",fontsize=12,fontweight="bold",y=0.985)
hs,ls=[],[]
for key,title,lab,psi in SC:
    d=RES[key]; xe=d["e"]/R
    h,=b0.step(np.r_[xe[:-1],xe[-1]],np.r_[d["rho"],d["rho"][-1]],where="post",color=COLS[key],lw=1.7); hs.append(h)
    ls.append(f"{title}  ({lab}):  densest {d['rho'].max():.3g}, emptiest {d['rho'].min():.3g}  ({d['rho'].max()/d['rho'].min():.0f}×)")
    b1.step(np.r_[xe[:-1],xe[-1]],np.r_[d["w"],d["w"][-1]],where="post",color=COLS[key],lw=1.7)
rr=M/vol(eref)/RHO_REF0
h,=b0.step(np.r_[eref[:-1],eref[-1]]/R,np.r_[rr,rr[-1]],where="post",color="0.4",lw=1.1,ls="--"); hs.append(h); ls.append(f"unmodulated (A = 0, width {eref[1]-eref[0]:.2f}):  densest 1, emptiest {rr[-1]:.3g}  ({1/rr[-1]:.0f}×)")
b1.axhline(eref[1]-eref[0],color="0.4",lw=1.1,ls="--")
b0.set_yscale("log"); b0.set_ylabel("shell density ÷ innermost shell\nof the unmodulated set"); b0.grid(alpha=0.3,which="both")
b1.set_ylabel("shell width  (radius units, R = 15)"); b1.grid(alpha=0.3); b1.set_xlim(0,1); b1.set_xlabel("normalised radius  r / R")
for a in (b0,b1): a.axvline(R0/R,color="0.6",lw=0.8)
b0.set_title("(a) shell density",loc="left",fontsize=10,fontweight="bold"); b1.set_title("(b) shell width",loc="left",fontsize=10,fontweight="bold")
fig.legend(hs,ls,loc="upper center",bbox_to_anchor=(0.53,0.135),ncol=1,frameon=False,fontsize=9)
fig.savefig("resonant_compare.png",dpi=250); plt.close(fig)
info={k:dict(c=v["c"],widths=v["w"].round(4).tolist(),width_min=float(v["w"].min()),width_max=float(v["w"].max()),rho=v["rho"].round(5).tolist(),
      rho_max=float(v["rho"].max()),rho_min=float(v["rho"].min()),ratio=float(v["rho"].max()/v["rho"].min()),fit_max_rel_err=v["fit_max_rel_err"]) for k,v in RES.items()}
info["_meta"]=dict(N=N,r_in=R0,R=R,A=A,norm=[NORM.vmin,NORM.vmax],ref_outer=float(rr[-1]),page_in=[W,H])
json.dump(info,open("resonant.json","w"),indent=1)
for k,v in info.items(): print(k,{a:b for a,b in v.items() if a not in("widths","rho")})
