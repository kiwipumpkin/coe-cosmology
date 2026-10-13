"""Pages 3 and 4 (Chris, 11 Oct): same three spacings and radii as shells_remake.py / structure_growth.py, but UNEQUAL mass per shell,
mass decreasing inward (outermost shell holds the most):
    M_i = M_tot * q_m^i / sum_j q_m^j ,  i = 0 (innermost) .. 11,   q_m = MASS_RATIO^(1/(N-1))   (MASS_RATIO = M_11/M_0 = 3 -> q_m = 1.1050)
    M_tot = 12 M (same total as the equal-mass case).
    rho_i = M_i / [(4/3) pi (r_{i+1}^3 - r_i^3)], shown DIVIDED BY THE INNERMOST-SHELL DENSITY OF THE EQUAL-MASS CASE with the same spacing
    (so colours compare directly with pages 1-2; the inner shell here is 0.544).
Smooth fit: rho_fit(x) = A x^-n + C x^p   (x = r/R; C >= 0; p free). p = 0 is the old 'A x^-n + C' (flattening); p > 0 lets it RISE again
    (needed for narrowing, whose density has a minimum at shell 8 and then rises). Least squares in log rho at shell midpoints.
Wave: cos(phi), phi = 2 pi x (shells crossed) -- unchanged (depends only on spacing). Bottom: rho_fit(x)/rho_fit(x_in) * cos(phi), +/- envelope.
Structure: same rules as structure_growth.py:
    sigma_i = rms(rho/rho_i - 1) = 0.8 * i/(N-1) * (rho_i / rho_i,const-equal-mass)^(1/3), first shell smooth, monotonic rise (asserted);
    spacing lambda_i = 0.6 w_i; lognormal field with mean over the shell's drawn area exactly rho_i; same noise seed as pages 1-2.
Colour scales: identical to the existing pages (waves pages LogNorm 10^-2.1..1, structure pages LogNorm 10^-2.6..10^0.4); not extended
    (all smooth densities fit; only a few clump peaks in the narrowing outer shells saturate, fraction reported).
"""
import json,sys,numpy as np,matplotlib,matplotlib.ticker
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.colors import LogNorm,LinearSegmentedColormap
from scipy.optimize import brentq,curve_fit
from PIL import Image
plt.rcParams.update({"font.size":12})
N=12; R0=3.0; W0=1.0; SPAN=45.0; MASS_RATIO=3.0; S_OUT=0.8; LFAC=0.6; BW=0.35; DPI=130
QM=MASS_RATIO**(1/(N-1))
CMAP=LinearSegmentedColormap.from_list("shells",["#08154f","#1747a6","#4f86cf","#8a6db0","#b45a8e","#d45a6a","#9b0f1a"])
NORM_W=LogNorm(10**-2.1,1.0); NORM_S=LogNorm(10**-2.6,10**0.4)
CASES=[("constant","Constant shell spacing",15.00),("widening","Widening shell spacing (geometric growth)",24.00),("narrowing","Narrowing shell spacing (geometric shrink)",12.37)]
def edges(q): return R0+np.r_[0,np.cumsum(W0*q**np.arange(N))]
def q_for(Rout): return 1.0 if abs(Rout-15)<1e-9 else brentq(lambda q: edges(q)[-1]-Rout,0.5,2.0)
Mi=QM**np.arange(N); Mi=N*Mi/Mi.sum()
def dens(e,masses):
    V=4/3*np.pi*(e[1:]**3-e[:-1]**3); return masses/V/(1/V[0])          # / equal-mass inner density
RHOC=dens(edges(1.0),np.ones(N))
def fitrho(xe,rho):
    xm=0.5*(xe[1:]+xe[:-1]); f=lambda x,lA,n,lC,p: np.log(np.exp(lA)*x**(-n)+np.exp(lC)*x**p)
    best=None
    for p0 in (0.0,1.0,3.0):
        try:
            P,_=curve_fit(f,xm,np.log(rho),p0=[np.log(rho[0]*xm[0]**2),2.0,np.log(rho.min()),p0],bounds=([-30,0,-30,-1],[30,12,10,12]),maxfev=40000)
            r=np.sum((f(xm,*P)-np.log(rho))**2)
            if best is None or r<best[0]: best=(r,P)
        except RuntimeError: pass
    P=best[1]; A,n,C,p=np.exp(P[0]),P[1],np.exp(P[2]),P[3]
    return (lambda x: A*x**(-n)+C*x**p),dict(A=float(A),n=float(n),C=float(C),p=float(p),max_rel_err=float(np.exp(np.abs(f(xm,*P)-np.log(rho))).max()-1))
def fitlabel(pf): return rf"smooth fit $\rho={pf['A']:.3g}\,x^{{-{pf['n']:.2f}}}+{pf['C']:.3g}\,x^{{{pf['p']:.2f}}}$,  x = r/R"
def waves(e,x):
    r=x*e[-1]; i=np.clip(np.searchsorted(e,r,side="right")-1,0,N-1); phi=2*np.pi*(i+(r-e[i])/(e[i+1]-e[i])); return np.where(r>=e[0],np.cos(phi),np.nan)
TXT=lambda Rout,q: f"{N} shells, mass rises outward ×{QM:.4f} per shell\n(outermost = {MASS_RATIO:g}× innermost; same total mass)\ninner radius {R0:g}, outer radius R = {Rout:.2f}"+("" if q==1 else f"\neach shell {q:.4f}× as wide as the one inside it")
# ---------------- layout helper ----------------
LM=1.1; AW=8.6; Wf=10.0; HW=AW*np.sin(np.radians(SPAN))*1.01
def frame(title,hs,cbar):
    hcb,gcb=(0.22,1.05) if cbar else (0,0.75); base=1.0 if cbar else 0.0; top=0.95
    H=base+hcb+gcb+sum(hs)+0.22*(len(hs)-1)+HW+top; fig=plt.figure(figsize=(Wf,H)); ax_=lambda y,h: fig.add_axes([LM/Wf,y/H,AW/Wf,h/H])
    fig.suptitle(title,fontsize=15,fontweight="bold",y=1-0.3/H)
    y=base+hcb+gcb; axs=[]
    for h in hs[::-1]: axs.append(ax_(y,h)); y+=h+0.22
    y-=0.22; aw=ax_(y,HW); aw.set_xlim(0,1); aw.set_ylim(0,HW/AW); aw.axis("off")
    cax=fig.add_axes([LM/Wf,base/H,AW/Wf,hcb/H]) if cbar else None
    return fig,aw,axs[::-1],cax
def finish(axs,xe,cax,norm,label):
    for a in axs:
        for b in xe: a.axvline(b,color="0.8",lw=0.8,zorder=0)
        a.set_xlim(0,1); a.grid(axis="y",alpha=0.25)
    for a in axs[:-1]: plt.setp(a.get_xticklabels(),visible=False)
    axs[-1].set_xlabel("normalised radius  r / R   (grey lines = shell boundaries)")
    if cax is not None:
        cb=plt.gcf().colorbar(plt.cm.ScalarMappable(norm=norm,cmap=CMAP),cax=cax,orientation="horizontal"); cb.set_label(label,fontsize=11)
CBL_W="shell density ÷ innermost shell of the equal-mass case (same logarithmic colour scale as pages 1–2)"
CBL_S="density ÷ innermost shell of the equal-mass case (one logarithmic colour scale, same as page 2)"
# ---------------- waves (page 3) ----------------
def wave_fig(key,title,Rout,cbar,out):
    q=q_for(Rout); e=edges(q); xe=e/e[-1]; rho=dens(e,Mi); g,pf=fitrho(xe,rho)
    fig,aw,(a1,a2,a3),cax=frame(f"Outer-heavy concentric shells — {title}",[2.5,2.0,2.0],cbar)
    th=np.linspace(0,np.radians(SPAN),300); band=lambda r0,r1: np.r_[np.c_[r1*np.cos(th),r1*np.sin(th)],np.c_[r0*np.cos(th[::-1]),r0*np.sin(th[::-1])]]
    aw.add_patch(Polygon(np.r_[[[0,0]],np.c_[xe[0]*np.cos(th),xe[0]*np.sin(th)]],fc="#dcdcdc",ec="none"))
    for k in range(N): aw.add_patch(Polygon(band(xe[k],xe[k+1]),fc=CMAP(NORM_W(rho[k])),ec="white",lw=1.0))
    aw.text(0.0,HW/AW*0.97,TXT(Rout,q),fontsize=12,va="top",ha="left")
    x=np.linspace(0,1,4000); wv=waves(e,x); rf=np.where(x>=xe[0],g(np.maximum(x,1e-9)),np.nan); env=rf/g(xe[0])
    rhoeq=dens(e,np.ones(N))
    a1.step(np.r_[xe[:-1],xe[-1]],np.r_[rhoeq,rhoeq[-1]],where="post",color="0.6",lw=1.0,ls=":",label="equal-mass case (pages 1–2)")
    for k in range(N): a1.hlines(rho[k],xe[k],xe[k+1],color=CMAP(NORM_W(rho[k])),lw=4)
    a1.plot(x,rf,"--",color="0.15",lw=1.3,label=fitlabel(pf))
    a1.set_yscale("log"); a1.set_ylim(10**-2.3,1.6); a1.set_ylabel("shell density\n(÷ equal-mass inner)"); a1.legend(loc="upper right",fontsize=9.5,frameon=False)
    a1.text(0.015,0.06,f"physical outer radius R = {Rout:.2f}",transform=a1.transAxes,fontsize=11,color="0.3")
    a2.plot(x,wv,color="#2c3e50",lw=1.6); a2.set_ylabel("wave following\nshell spacing"); a2.set_ylim(-1.25,1.25)
    a3.plot(x,env,color="#b03050",lw=1.0,ls="--",label="± smoothed density fit"); a3.plot(x,-env,color="#b03050",lw=1.0,ls="--")
    a3.plot(x,env*wv,color="#2c3e50",lw=1.6,label="smoothed density fit × wave"); a3.set_ylabel("smoothed density\nfit × wave"); a3.set_ylim(-1.25,1.25)
    a3.legend(loc="upper right",fontsize=10,frameon=False,ncol=2)
    for a in (a2,a3): a.axhline(0,color="0.6",lw=0.8)
    finish([a1,a2,a3],xe,cax,NORM_W,CBL_W); fig.savefig(out,dpi=DPI); plt.close(fig)
    return dict(q=q,mass=Mi.tolist(),rho=rho.tolist(),rho_equal_mass=rhoeq.tolist(),fit=pf)
# ---------------- structure (page 4) ----------------
def sbuild(e,rho,nx=1100,seed=11):
    w=np.diff(e); Rout=e[-1]; sig=S_OUT*np.arange(N)/(N-1)*(rho/RHOC)**(1/3); lam=LFAC*w
    assert np.all(np.diff(sig)>0),sig
    ny=int(nx*np.sin(np.radians(SPAN)))+2; dx=Rout/nx; X=(np.arange(nx)+0.5)*dx; Y=(np.arange(ny)+0.5)*dx
    XX,YY=np.meshgrid(X,Y); RR=np.hypot(XX,YY); TT=np.degrees(np.arctan2(YY,XX))
    px,py=nx+512,ny+512; K=np.hypot(*np.meshgrid(2*np.pi*np.fft.fftfreq(px,dx),2*np.pi*np.fft.fftfreq(py,dx)))
    Wn=np.fft.fft2(np.random.default_rng(seed).normal(size=(py,px))); D=np.full((ny,nx),np.nan); rows=[]
    kb=np.geomspace(2*np.pi/(px*dx)*2,np.pi/dx,240); idx=np.digitize(K.ravel(),kb); cnt=np.maximum(np.bincount(idx,minlength=len(kb)+1)[1:len(kb)],1)
    for i in range(N):
        m=(RR>=e[i])&(RR<e[i+1])&(TT<=SPAN)
        if sig[i]==0: D[m]=rho[i]; rows.append(dict(shell=1,rho=float(rho[i]),sigma_target=0.,sigma_measured=0.,spacing_target=float(lam[i]),spacing_measured=None,mean_over_rho=1.0,max_over_rho=1.0,min_over_rho=1.0,pct={k_:float(rho[i]) for k_ in ("min","p2.5","p16","p84","p97.5","max")})); continue
        ki=2*np.pi/lam[i]; T=np.where(K>0,np.exp(-np.log(np.maximum(K,1e-12)/ki)**2/(4*BW**2))/np.maximum(K,1e-12),0.)
        G=np.real(np.fft.ifft2(Wn*T)); g=G[:ny,:nx]; g=(g-g[m].mean())/g[m].std()
        s=np.sqrt(np.log(1+sig[i]**2)); L=np.exp(s*g-s*s/2); L=L/L[m].mean(); D[m]=rho[i]*L[m]
        pw=np.bincount(idx,weights=(np.abs(np.fft.fft2(G))**2*K**2).ravel(),minlength=len(kb)+1)[1:len(kb)]/cnt
        kpk=np.sqrt(kb[:-1]*kb[1:])[np.argmax(pw)]; r_=L[m]-1
        rows.append(dict(shell=i+1,rho=float(rho[i]),sigma_target=float(sig[i]),sigma_measured=float(np.sqrt(np.mean(r_**2))),spacing_target=float(lam[i]),
            spacing_measured=float(2*np.pi/kpk),mean_over_rho=float(L[m].mean()),max_over_rho=float(L[m].max()),min_over_rho=float(L[m].min()),pct=dict(zip(("min","p2.5","p16","p84","p97.5","max"),[float(v) for v in rho[i]*np.r_[L[m].min(),np.percentile(L[m],[2.5,16,84,97.5]),L[m].max()]]))))
    return dict(rows=rows,D=D,core=(RR<e[0])&(TT<=SPAN),asp=ny/nx)
def struct_fig(key,title,Rout,cbar,out,B=None):
    q=q_for(Rout); e=edges(q); xe=e/e[-1]; rho=dens(e,Mi); g,pf=fitrho(xe,rho)
    if B is None: B=sbuild(e,rho)
    rows=B["rows"]; fig,aw,(a1,a2,a3),cax=frame(f"Structure in outer-heavy shells — {title}",[2.4,2.1,2.1],cbar)
    img=CMAP(NORM_S(np.clip(np.nan_to_num(B["D"],nan=1.0),NORM_S.vmin,NORM_S.vmax))); img[np.isnan(B["D"])]=(1,1,1,0); img[B["core"]]=(0.86,0.86,0.86,1)
    aw.imshow(img,origin="lower",extent=[0,1,0,B["asp"]],interpolation="nearest",zorder=1)
    th=np.linspace(0,np.radians(SPAN),300)
    for b in xe: aw.plot(b*np.cos(th),b*np.sin(th),color="white",lw=0.9,zorder=3)
    aw.text(0.0,HW/AW*0.97,f"mass rises outward (outermost = {MASS_RATIO:g}× innermost)\ninnermost shell smooth; clumpier outward\neach shell's clumpy density averages to its shell density",fontsize=11.5,va="top",ha="left",zorder=4)
    x=np.linspace(xe[0],1,800)
    from matplotlib.patches import Patch
    for k in range(N):
        p=rows[k]["pct"]; c=CMAP(NORM_S(rho[k]))
        a1.fill_between([xe[k],xe[k+1]],p["p2.5"],p["p97.5"],color=c,alpha=0.16,lw=0,zorder=1)
        a1.fill_between([xe[k],xe[k+1]],p["p16"],p["p84"],color=c,alpha=0.42,lw=0,zorder=2)
        if k>0: a1.errorbar([0.5*(xe[k]+xe[k+1])],[rho[k]],yerr=[[rho[k]-p["min"]],[p["max"]-rho[k]]],fmt="none",ecolor="0.25",elinewidth=0.9,capsize=3.5,capthick=0.9,zorder=3)
    for k in range(N): a1.hlines(rho[k],xe[k],xe[k+1],color=CMAP(NORM_S(rho[k])),lw=4,zorder=4)
    a1.plot(x,g(x),"--",color="0.15",lw=1.2,label=fitlabel(pf))
    xm=0.5*(xe[1:]+xe[:-1]); a1.plot(xm,[r["mean_over_rho"]*r["rho"] for r in rows],"o",zorder=5,ms=5,mfc="white",mec="k",mew=1.0,label="mean of the clumpy field in each shell (= shell density)")
    a1.set_yscale("log"); a1.set_ylim(10**-3.2,200.0); a1.set_ylabel("shell density\n(÷ equal-mass inner)")
    h_,l_=a1.get_legend_handles_labels(); from matplotlib.lines import Line2D
    h_+=[Patch(color="0.45",alpha=0.42),Patch(color="0.45",alpha=0.16),Line2D([],[],color="0.25",lw=0.9,marker="_",ms=7)]; l_+=["clump spread: 16–84% of density in the shell","clump spread: 2.5–97.5%","max–min density in the shell"]
    a1.legend(h_,l_,loc="upper right",fontsize=9,frameon=False,ncol=2)
    sm=[r["sigma_measured"] for r in rows]; sp=[r["spacing_measured"] for r in rows]
    for k in range(N): a2.hlines(sm[k],xe[k],xe[k+1],color="#7b241c",lw=3.5)
    a2.set_ylim(0,1.45); a2.set_ylabel("clumpiness\nrms(ρ/ρ_shell − 1)")
    for k in range(1,N): a3.hlines(sp[k],xe[k],xe[k+1],color="#1f3a93",lw=3.5)
    a3.text(0.5*(xe[0]+xe[1]),0.3,"smooth",ha="center",fontsize=9,color="0.4")
    a3.set_yscale("log"); a3.set_ylim(0.2,3); a3.set_yticks([0.3,0.5,1,2]); a3.set_yticklabels(["0.3","0.5","1","2"]); a3.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a3.set_ylabel("typical clump spacing\n(radius units)")
    finish([a1,a2,a3],xe,cax,NORM_S,CBL_S); fig.savefig(out,dpi=DPI); plt.close(fig)
    sat=float(np.mean(B["D"][~np.isnan(B["D"])]>NORM_S.vmax)); low=float(np.mean(B["D"][~np.isnan(B["D"])]<NORM_S.vmin))
    return B,dict(fit=pf,rows=rows,frac_pixels_above_scale=sat,frac_pixels_below_scale=low)
def cbar_strip(norm,label,out):
    fig=plt.figure(figsize=(Wf,1.5)); cax=fig.add_axes([LM/Wf,0.62,AW/Wf,0.22/1.5])
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=CMAP),cax=cax,orientation="horizontal"); cb.set_label(label,fontsize=12); fig.savefig(out,dpi=DPI); plt.close(fig)
def stack(files,out):
    ims=[Image.open(f).convert("RGB") for f in files]; W=max(i.width for i in ims); H=sum(i.height for i in ims)
    P=Image.new("RGB",(W,H),"white"); y=0
    for i in ims: P.paste(i,(0,y)); y+=i.height
    P.save(out+".png"); P.save(out+".pdf",resolution=150.0); return P.size
if __name__=="__main__":
    import os; os.makedirs("_tmp_pages",exist_ok=True); info=dict(mass_ratio=MASS_RATIO,q_m=QM,masses=Mi.tolist(),waves={},structure={})
    for key,title,Rout in CASES:
        info["waves"][key]=wave_fig(key,title,Rout,True,f"shells_massout_{key}.png"); wave_fig(key,title,Rout,False,f"_tmp_pages/w_{key}.png")
        B,inf=struct_fig(key,title,Rout,True,f"structure_massout_{key}.png"); struct_fig(key,title,Rout,False,f"_tmp_pages/s_{key}.png",B); info["structure"][key]=inf
    cbar_strip(NORM_W,CBL_W,"_tmp_pages/cb_w.png"); cbar_strip(NORM_S,CBL_S,"_tmp_pages/cb_s.png")
    print(stack([f"_tmp_pages/w_{k}.png" for k,_,_ in CASES]+["_tmp_pages/cb_w.png"],"page3_shells_waves_massout"))
    print(stack([f"_tmp_pages/s_{k}.png" for k,_,_ in CASES]+["_tmp_pages/cb_s.png"],"page4_shells_structure_massout"))
    json.dump(info,open("shells_massout.json","w"),indent=1)
    for k in info["waves"]: print(k,"fit",{a:round(b,4) for a,b in info["waves"][k]["fit"].items()},"sat",info["structure"][k]["frac_pixels_above_scale"],info["structure"][k]["frac_pixels_below_scale"])
    for k in info["structure"]: print(k,[round(r["sigma_measured"],3) for r in info["structure"][k]["rows"]],[None if r["spacing_measured"] is None else round(r["spacing_measured"],2) for r in info["structure"][k]["rows"]])
