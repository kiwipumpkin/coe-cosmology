"""Remake of Chris's 'Equal-Mass Concentric Shells' figure (11 Oct 2026), one PNG per case + a density comparison.
Physics
  12 shells, inner radius r_0 = 3 (figure units), first shell width w_0 = 1, every shell holds the same mass M:
      rho_i = M / [ (4/3) pi (r_{i+1}^3 - r_i^3) ]            (shown normalised to the innermost shell, rho_i / rho_0)
  spacings: constant  w_i = w_0                          -> R = 3 + 12   = 15.00
            widening  w_i = w_0 q^i, q > 1 (geometric)   -> q solved so R = 24.00
            narrowing w_i = w_0 q^i, q < 1 (geometric)   -> q solved so R = 12.37
  (r_0 = 3 and w_0 = 1 are read off the original: r_in/R = 0.200, 0.125, 0.243 and first widths ~1 in all three.)
  wave following the shell spacing: local wavelength lambda(r) = width of the shell containing r,
      phi(r) = 2 pi * integral_{r_0}^{r} dr'/lambda(r')  (= 2 pi x number of shells crossed), wave = cos(phi)  -> crests on every boundary
  smoothed density fit x wave (Chris, 11 Oct; replaces the earlier j0 panel):
      rho_fit(x) = A x^(-n) + C,  x = r/R, fitted by least squares in log(rho) to the shell densities at the shell midpoints
      (rho_i/rho_0, all shells equally weighted; A, n, C >= 0). The constant C lets the profile flatten (narrowing case).
      bottom panel: rho_fit(x)/rho_fit(x_in) * cos(phi(x)), same phi as the wave panel, with +/- rho_fit/rho_fit(x_in) dashed.
"""
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from scipy.optimize import brentq, curve_fit
plt.rcParams.update({"font.size":12})
N=12; R0=3.0; W0=1.0; M=1.0
CMAP=LinearSegmentedColormap.from_list("shells",["#08154f","#1747a6","#4f86cf","#8a6db0","#b45a8e","#d45a6a","#9b0f1a"])
NORM=LogNorm(10**-2.1,1.0)
def edges(q):
    w=W0*q**np.arange(N); return R0+np.r_[0,np.cumsum(w)]
def q_for(Rout): return brentq(lambda q: edges(q)[-1]-Rout,0.5,2.0) if abs(Rout-(R0+N*W0))>1e-9 else 1.0
CASES=[("constant","Constant shell spacing",15.00),("widening","Widening shell spacing (geometric growth)",24.00),("narrowing","Narrowing shell spacing (geometric shrink)",12.37)]
def build(Rout):
    q=q_for(Rout); e=edges(q); rho=M/(4/3*np.pi*(e[1:]**3-e[:-1]**3)); return q,e,rho/rho[0]
def waves(e,x):
    r=x*e[-1]; i=np.clip(np.searchsorted(e,r,side="right")-1,0,N-1); frac=(r-e[i])/(e[i+1]-e[i])
    phi=2*np.pi*(i+frac); P0=2*np.pi*R0/W0+np.pi/2; Phi=phi+P0
    j0=np.sin(Phi)/Phi/(np.sin(P0)/P0); inside=r>=e[0]
    return np.where(inside,np.cos(phi),np.nan),np.where(inside,j0,np.nan)
def fitrho(xe,rho):
    xm=0.5*(xe[1:]+xe[:-1]); f=lambda x,lA,n,C: np.log(np.exp(lA)*x**(-n)+C)
    p,_=curve_fit(f,xm,np.log(rho),p0=[np.log(rho[0]*xm[0]**2),2.0,1e-3],bounds=([-30,0,0],[30,10,10]),maxfev=20000)
    A,n,C=np.exp(p[0]),p[1],p[2]; g=lambda x: A*x**(-n)+C
    err=np.exp(np.abs(f(xm,*p)-np.log(rho))).max()-1
    return g,dict(A=float(A),n=float(n),C=float(C),max_rel_err=float(err))
FITS={}
out={}
for key,title,Rout in CASES:
    q,e,rho=build(Rout); xe=e/e[-1]; out[key]=(q,e,rho)
    # explicit layout in inches so the wedge has EQUAL aspect and the same x-scale (r/R) as the panels below
    SPAN=45.0; W=10.0; LM=1.1; AW=8.6; HW=AW*np.sin(np.radians(SPAN))*1.01
    hcb,gcb,h3,g,h2,h1=0.22,1.05,2.0,0.22,2.0,2.5; base=0.8; top=0.95
    H=base+hcb+gcb+h3+g+h2+g+h1+HW+top
    fig=plt.figure(figsize=(W,H)); ax_=lambda y,h: fig.add_axes([LM/W,y/H,AW/W,h/H])
    y3=base+hcb+gcb; y2=y3+h3+g; y1=y2+h2+g; yw=y1+h1
    fig.suptitle(f"Equal-mass concentric shells — {title}",fontsize=15,fontweight="bold",y=1-0.3/H)
    aw=ax_(yw,HW); aw.set_xlim(0,1); aw.set_ylim(0,HW/AW); aw.axis("off")
    th=np.linspace(0,np.radians(SPAN),300)
    from matplotlib.patches import Polygon
    def band(r0,r1): return np.r_[np.c_[r1*np.cos(th),r1*np.sin(th)],np.c_[r0*np.cos(th[::-1]),r0*np.sin(th[::-1])]]
    aw.add_patch(Polygon(np.r_[[[0,0]],np.c_[xe[0]*np.cos(th),xe[0]*np.sin(th)]],fc="#dcdcdc",ec="none"))
    for k in range(N): aw.add_patch(Polygon(band(xe[k],xe[k+1]),fc=CMAP(NORM(rho[k])),ec="white",lw=1.0))
    aw.text(0.0,HW/AW*0.97,f"{N} shells, equal mass each\ninner radius {R0:g}, outer radius R = {Rout:.2f}"+("" if q==1 else f"\neach shell {q:.4f}× as wide\nas the one inside it"),fontsize=12,va="top",ha="left")
    a1=ax_(y1,h1); a2=ax_(y2,h2); a3=ax_(y3,h3); cax=fig.add_axes([LM/W,base/H,AW/W,hcb/H])
    x=np.linspace(0,1,4000); wv,_=waves(e,x); g,pf=fitrho(xe,rho); FITS[key]=pf
    xin=x>=xe[0]; rf=np.where(xin,g(np.maximum(x,1e-9)),np.nan); env=rf/g(xe[0]); dfit=env*wv
    for k in range(N): a1.hlines(rho[k],xe[k],xe[k+1],color=CMAP(NORM(rho[k])),lw=4)
    a1.set_yscale("log"); a1.set_ylim(10**-2.3,1.6); a1.set_ylabel("shell density\n(÷ innermost shell)")
    a1.text(0.015,0.06,f"physical outer radius R = {Rout:.2f}",transform=a1.transAxes,fontsize=11,color="0.3")
    a2.plot(x,wv,color="#2c3e50",lw=1.6); a2.set_ylabel("wave following\nshell spacing"); a2.set_ylim(-1.25,1.25)
    a3.plot(x,env,color="#b03050",lw=1.0,ls="--",label="± smoothed density fit"); a3.plot(x,-env,color="#b03050",lw=1.0,ls="--")
    a3.plot(x,dfit,color="#2c3e50",lw=1.6,label="smoothed density fit × wave"); a3.set_ylabel("smoothed density\nfit × wave"); a3.set_ylim(-1.25,1.25)
    a3.legend(loc="upper right",fontsize=10,frameon=False,ncol=2)
    a1.plot(x,rf,color="0.15",lw=1.3,ls="--",label=(rf"smooth fit  $\rho = {pf['A']:.3g}\,(r/R)^{{-{pf['n']:.2f}}}$"+(rf"$\,+\,{pf['C']:.3g}$" if pf['C']>1e-6 else "")))
    a1.legend(loc="upper right",fontsize=10,frameon=False)
    for a in (a1,a2,a3):
        for b in xe: a.axvline(b,color="0.8",lw=0.8,zorder=0)
        a.set_xlim(0,1); a.grid(axis="y",alpha=0.25)
    for a in (a2,a3): a.axhline(0,color="0.6",lw=0.8)
    plt.setp(a1.get_xticklabels(),visible=False); plt.setp(a2.get_xticklabels(),visible=False)
    a3.set_xlabel("normalised radius  r / R   (grey lines = shell boundaries)")
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=NORM,cmap=CMAP),cax=cax,orientation="horizontal")
    cb.set_label("shell density ÷ innermost shell (same logarithmic colour scale in all three figures)")
    fig.savefig(f"shells_{key}.png",dpi=130); plt.close(fig)
# comparison
fig,ax=plt.subplots(figsize=(10,6.5)); cols={"constant":"#555555","widening":"#1747a6","narrowing":"#b03050"}
for key,title,Rout in CASES:
    q,e,rho=out[key]; xe=e/e[-1]
    ax.step(np.r_[xe[:-1],xe[-1]],np.r_[rho,rho[-1]],where="post",color=cols[key],lw=2.4,label=f"{key}: R = {Rout:.2f}"+("" if q==1 else f", width ratio {q:.4f}")+f", outer/inner density {rho[-1]:.3f}")
ax.set_yscale("log"); ax.set_xlim(0,1); ax.set_xlabel("normalised radius  r / R"); ax.set_ylabel("shell density ÷ innermost shell"); ax.grid(alpha=0.3,which="both")
ax.set_title("Equal-mass shells: density profile for the three spacings",fontsize=14,fontweight="bold")
ax.legend(loc="upper center",bbox_to_anchor=(0.5,-0.14),ncol=1,frameon=False,fontsize=11)
fig.savefig("shells_compare_density.png",dpi=130,bbox_inches="tight"); plt.close(fig)
import json; json.dump(FITS,open("shells_fits.json","w"),indent=1); print(FITS)
for k,(q,e,rho) in out.items(): print(k,"q=%.5f"%q,"edges/R",np.round(e/e[-1],3).tolist(),"rho",np.round(rho,4).tolist())
