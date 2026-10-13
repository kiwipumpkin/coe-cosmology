"""Structure in the three equal-mass shell cases, v2 (Chris, 11 Oct). Shells = shells_remake.py (12 shells, inner radius 3, first width 1;
constant / widening q=1.0970 / narrowing q=0.9537; R = 15 / 24 / 12.37; equal mass per shell).
RULES (Chris):
 colour      : ONE log colour scale of actual density relative to the innermost shell (same gradient as the shell images), identical in all three.
 clumpiness  : shell 1 (innermost) perfectly smooth; clumpiness rises monotonically outward (outer shells seen later / developed longer).
               sigma_i = rms(rho/rho_i - 1) = S_OUT * (i/(N-1)) * (rho_i / rho_i,const)^(1/3),   i = 0..11
               S_OUT = 0.8 (outermost shell of the constant case); rho_i,const = density of the same shell in the constant case,
               so constant = the base ramp, widening (diluted) is fainter, narrowing (concentrated) is sharper. Exponent 1/3 keeps the
               rise monotonic in every case (checked).
 clump size  : typical clump spacing lambda_i = 0.6 x shell width w_i  (widening -> larger clumps, narrowing -> smaller, constant -> fixed).
 field       : Gaussian random field g with power concentrated around k_i = 2 pi / lambda_i (log-normal band, width 0.35 in ln k), same white-noise
               seed for every shell; lognormal density rho = rho_i * exp(s g - s^2/2), s = sqrt(ln(1 + sigma_i^2)), then rescaled so that the
               mean over the shell's drawn area is EXACTLY rho_i (mass conservation; positive everywhere).
 measured    : clumpiness = rms(rho/rho_i - 1) over the shell's pixels; spacing = 2 pi / k at the peak of the field's power spectrum.
"""
import json,numpy as np,matplotlib,matplotlib.ticker
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm,LinearSegmentedColormap
from scipy.optimize import brentq,curve_fit
plt.rcParams.update({"font.size":12})
N=12; R0=3.0; W0=1.0; S_OUT=0.8; LFAC=0.6; BW=0.35; SPAN=45.0
CASES=[("constant","Constant shell spacing",15.00),("widening","Widening shell spacing (geometric growth)",24.00),("narrowing","Narrowing shell spacing (geometric shrink)",12.37)]
CMAP=LinearSegmentedColormap.from_list("shells",["#08154f","#1747a6","#4f86cf","#8a6db0","#b45a8e","#d45a6a","#9b0f1a"])
NORM=LogNorm(10**-2.6,10**0.4)     # same gradient as the shell images; range widened from 10^-2.1..1 to hold clump peaks and troughs
def edges(q): return R0+np.r_[0,np.cumsum(W0*q**np.arange(N))]
def q_for(Rout): return 1.0 if abs(Rout-15)<1e-9 else brentq(lambda q: edges(q)[-1]-Rout,0.5,2.0)
def rho_of(e): r=1/(4/3*np.pi*(e[1:]**3-e[:-1]**3)); return r/r[0]
def fitrho(xe,rho):
    xm=0.5*(xe[1:]+xe[:-1]); f=lambda x,lA,n,C: np.log(np.exp(lA)*x**(-n)+C)
    p,_=curve_fit(f,xm,np.log(rho),p0=[np.log(rho[0]*xm[0]**2),2.0,1e-3],bounds=([-30,0,0],[30,10,10]),maxfev=20000)
    return (lambda x: np.exp(p[0])*x**(-p[1])+p[2]),dict(A=float(np.exp(p[0])),n=float(p[1]),C=float(p[2]))
RHOC=rho_of(edges(1.0))
def build(key,Rout,nx=1100,seed=11):
    q=q_for(Rout); e=edges(q); rho=rho_of(e); w=np.diff(e)
    sig=S_OUT*np.arange(N)/(N-1)*(rho/RHOC)**(1/3); lam=LFAC*w
    assert np.all(np.diff(sig)>0), (key,sig)
    ny=int(nx*np.sin(np.radians(SPAN)))+2; dx=Rout/nx
    X=(np.arange(nx)+0.5)*dx; Y=(np.arange(ny)+0.5)*dx; XX,YY=np.meshgrid(X,Y); RR=np.hypot(XX,YY); TT=np.degrees(np.arctan2(YY,XX))
    # padded FFT grid (avoid wrap-around for the biggest clumps)
    px,py=nx+512,ny+512; kx=2*np.pi*np.fft.fftfreq(px,dx); ky=2*np.pi*np.fft.fftfreq(py,dx); K=np.hypot(*np.meshgrid(kx,ky))
    Wn=np.fft.fft2(np.random.default_rng(seed).normal(size=(py,px)))
    D=np.full((ny,nx),np.nan); rows=[]
    for i in range(N):
        m=(RR>=e[i])&(RR<e[i+1])&(TT<=SPAN)
        if sig[i]==0: D[m]=rho[i]; rows.append(dict(shell=1,rho=float(rho[i]),sigma_target=0.0,sigma_measured=0.0,spacing_target=float(lam[i]),spacing_measured=None,mean_over_rho=1.0,max_over_rho=1.0,min_over_rho=1.0,pct={k_:float(rho[i]) for k_ in ("min","p2.5","p16","p84","p97.5","max")})); continue
        ki=2*np.pi/lam[i]; T=np.where(K>0,np.exp(-np.log(np.maximum(K,1e-12)/ki)**2/(4*BW**2))/np.maximum(K,1e-12),0.0)   # power per log k ~ Gaussian in ln k (2D: P ~ T^2)
        G=np.real(np.fft.ifft2(Wn*T)); g=G[:ny,:nx]; g=(g-g[m].mean())/g[m].std()
        s=np.sqrt(np.log(1+sig[i]**2)); L=np.exp(s*g-s*s/2); L=L/L[m].mean()      # mean over the shell's drawn area = 1 exactly
        D[m]=rho[i]*L[m]
        # measured spacing: peak of the radially averaged power spectrum of g (full grid)
        Pk=np.abs(np.fft.fft2(G))**2; kb=np.geomspace(2*np.pi/(px*dx)*2,np.pi/dx,240); idx=np.digitize(K.ravel(),kb)
        pw=np.bincount(idx,weights=(Pk*K**2).ravel(),minlength=len(kb)+1)[1:len(kb)]/np.maximum(np.bincount(idx,minlength=len(kb)+1)[1:len(kb)],1)
        kpk=np.sqrt(kb[:-1]*kb[1:])[np.argmax(pw)]
        r_=L[m]-1; rows.append(dict(shell=i+1,rho=float(rho[i]),sigma_target=float(sig[i]),sigma_measured=float(np.sqrt(np.mean(r_**2))),spacing_target=float(lam[i]),
            spacing_measured=float(2*np.pi/kpk),mean_over_rho=float(L[m].mean()),max_over_rho=float(L[m].max()),min_over_rho=float(L[m].min()),pct=dict(zip(("min","p2.5","p16","p84","p97.5","max"),[float(v) for v in rho[i]*np.r_[L[m].min(),np.percentile(L[m],[2.5,16,84,97.5]),L[m].max()]]))))
    core=(RR<e[0])&(TT<=SPAN)
    return dict(q=q,e=e,rho=rho,rows=rows,D=D,core=core,asp=ny/nx)
RES={k:build(k,R) for k,_,R in CASES}
json.dump(dict(rules=__doc__,cases={k:dict(q=v["q"],R=float(v["e"][-1]),rows=v["rows"]) for k,v in RES.items()}),open("structure_growth.json","w"),indent=1)
def draw(key,title,Rout,cbar,out):
    B=RES[key]; e=B["e"]; xe=e/e[-1]; rho=B["rho"]; rows=B["rows"]; g,pf=fitrho(xe,rho)
    Wf=10.0; LM=1.1; AW=8.6; HW=AW*np.sin(np.radians(SPAN))*1.01
    hcb,gcb,h3,gp,h2,h1=(0.22 if cbar else 0),(1.05 if cbar else 0.75),2.1,0.22,2.1,2.4; base=(1.0 if cbar else 0.0); top=0.95; H=base+hcb+gcb+h3+gp+h2+gp+h1+HW+top
    fig=plt.figure(figsize=(Wf,H)); ax_=lambda y,h: fig.add_axes([LM/Wf,y/H,AW/Wf,h/H])
    y3=base+hcb+gcb; y2=y3+h3+gp; y1=y2+h2+gp; yw=y1+h1
    fig.suptitle(f"Structure in equal-mass shells — {title}",fontsize=15,fontweight="bold",y=1-0.3/H)
    aw=ax_(yw,HW); aw.set_xlim(0,1); aw.set_ylim(0,HW/AW); aw.axis("off")
    img=CMAP(NORM(np.clip(np.nan_to_num(B["D"],nan=1.0),NORM.vmin,NORM.vmax))); img[np.isnan(B["D"])]=(1,1,1,0); img[B["core"]]=(0.86,0.86,0.86,1)
    aw.imshow(img,origin="lower",extent=[0,1,0,B["asp"]],interpolation="nearest",zorder=1)
    th=np.linspace(0,np.radians(SPAN),300)
    for b in xe: aw.plot(b*np.cos(th),b*np.sin(th),color="white",lw=0.9,zorder=3)
    aw.text(0.0,HW/AW*0.97,"innermost shell smooth; clumpier outward\n(outer shells developed longer)\neach shell's clumpy density averages to its shell density",fontsize=11.5,va="top",ha="left",zorder=4)
    a1=ax_(y1,h1); a2=ax_(y2,h2); a3=ax_(y3,h3); cax=fig.add_axes([LM/Wf,base/H,AW/Wf,hcb/H]) if cbar else None
    x=np.linspace(xe[0],1,800)
    from matplotlib.patches import Patch
    for k in range(N):
        p=rows[k]["pct"]; c=CMAP(NORM(rho[k]))
        a1.fill_between([xe[k],xe[k+1]],p["p2.5"],p["p97.5"],color=c,alpha=0.16,lw=0,zorder=1)
        a1.fill_between([xe[k],xe[k+1]],p["p16"],p["p84"],color=c,alpha=0.42,lw=0,zorder=2)
        if k>0: a1.errorbar([0.5*(xe[k]+xe[k+1])],[rho[k]],yerr=[[rho[k]-p["min"]],[p["max"]-rho[k]]],fmt="none",ecolor="0.25",elinewidth=0.9,capsize=3.5,capthick=0.9,zorder=3)
    for k in range(N): a1.hlines(rho[k],xe[k],xe[k+1],color=CMAP(NORM(rho[k])),lw=4,zorder=4)
    a1.plot(x,g(x),"--",color="0.15",lw=1.2,label=rf"smooth fit $\rho={pf['A']:.3g}\,(r/R)^{{-{pf['n']:.2f}}}$"+(rf"$\,+\,{pf['C']:.3g}$" if pf['C']>1e-6 else ""))
    xm=0.5*(xe[1:]+xe[:-1]); mm=np.array([r["mean_over_rho"]*r["rho"] for r in rows])
    a1.plot(xm,mm,"o",zorder=5,ms=5,mfc="white",mec="k",mew=1.0,label="mean of the clumpy field in each shell (= shell density)")
    a1.set_yscale("log"); a1.set_ylim(10**-3.2,200.0); a1.set_ylabel("shell density\n(÷ innermost shell)")
    h_,l_=a1.get_legend_handles_labels(); from matplotlib.lines import Line2D
    h_+=[Patch(color="0.45",alpha=0.42),Patch(color="0.45",alpha=0.16),Line2D([],[],color="0.25",lw=0.9,marker="_",ms=7)]; l_+=["clump spread: 16–84% of density in the shell","clump spread: 2.5–97.5%","max–min density in the shell"]
    a1.legend(h_,l_,loc="upper right",fontsize=9,frameon=False,ncol=2)
    sm=np.array([r["sigma_measured"] for r in rows]); sp=[r["spacing_measured"] for r in rows]
    for k in range(N): a2.hlines(sm[k],xe[k],xe[k+1],color="#7b241c",lw=3.5)
    a2.set_ylim(0,1.45); a2.set_ylabel("clumpiness\nrms(ρ/ρ_shell − 1)")
    for k in range(1,N): a3.hlines(sp[k],xe[k],xe[k+1],color="#1f3a93",lw=3.5)
    a3.text(0.5*(xe[0]+xe[1]),0.3,"smooth",ha="center",fontsize=9,color="0.4")
    a3.set_yscale("log"); a3.set_ylim(0.2,3); a3.set_yticks([0.3,0.5,1,2]); a3.set_yticklabels(["0.3","0.5","1","2"]); a3.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    a3.set_ylabel("typical clump spacing\n(radius units)")
    for a_ in (a1,a2,a3):
        for b in xe: a_.axvline(b,color="0.8",lw=0.8,zorder=0)
        a_.set_xlim(0,1); a_.grid(axis="y",alpha=0.25)
    plt.setp(a1.get_xticklabels(),visible=False); plt.setp(a2.get_xticklabels(),visible=False)
    a3.set_xlabel("normalised radius  r / R   (grey lines = shell boundaries)")
    if cbar:
        cb=fig.colorbar(plt.cm.ScalarMappable(norm=NORM,cmap=CMAP),cax=cax,orientation="horizontal")
        cb.set_label(CBL,fontsize=11)
    fig.savefig(out,dpi=130); plt.close(fig)
CBL="density ÷ innermost shell (one logarithmic colour scale, identical in all three figures)"
import os,shells_massout as SM; os.makedirs("_tmp_p2",exist_ok=True)
for key,title,Rout in CASES: draw(key,title,Rout,True,f"structure_{key}.png"); draw(key,title,Rout,False,f"_tmp_p2/s_{key}.png")
SM.cbar_strip(NORM,"density ÷ innermost shell (one logarithmic colour scale, same for all three cases)","_tmp_p2/cb.png")
print("page2",SM.stack([f"_tmp_p2/s_{k}.png" for k,_,_ in CASES]+["_tmp_p2/cb.png"],"page2_shells_structure"))
import shutil; shutil.rmtree("_tmp_p2")
fig,ax=plt.subplots(2,1,figsize=(10,9),sharex=True); fig.subplots_adjust(hspace=0.12); cols={"constant":"#555555","widening":"#1747a6","narrowing":"#b03050"}
for key,title,Rout in CASES:
    B=RES[key]; xe=B["e"]/B["e"][-1]; s=[r["sigma_measured"] for r in B["rows"]]; sp=[r["spacing_measured"] or np.nan for r in B["rows"]]
    ax[0].step(np.r_[xe[:-1],xe[-1]],np.r_[s,s[-1]],where="post",color=cols[key],lw=2.4,label=f"{key} (R = {Rout:.2f})")
    ax[1].step(np.r_[xe[:-1],xe[-1]],np.r_[sp,sp[-1]],where="post",color=cols[key],lw=2.4)
ax[0].set_ylabel("clumpiness rms(ρ/ρ_shell − 1)"); ax[0].set_ylim(0,1.45); ax[0].grid(alpha=0.3)
ax[1].set_ylabel("typical clump spacing (radius units)"); ax[1].set_yscale("log"); ax[1].set_yticks([0.3,0.5,1,2]); ax[1].set_yticklabels(["0.3","0.5","1","2"]); ax[1].yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter()); ax[1].set_ylim(0.3,2.3); ax[1].grid(alpha=0.3,which="both"); ax[1].set_xlim(0,1); ax[1].set_xlabel("normalised radius  r / R")
ax[0].set_title("Structure in equal-mass shells: clumpiness and clump spacing vs radius",fontsize=14,fontweight="bold")
fig.legend(*ax[0].get_legend_handles_labels(),loc="lower center",ncol=3,frameon=False,fontsize=11,bbox_to_anchor=(0.5,0.0)); fig.subplots_adjust(bottom=0.13)
fig.savefig("structure_compare.png",dpi=130,bbox_inches="tight"); plt.close(fig)
for k,B in RES.items():
    print(k); [print("  %2d rho %.4f sig %.3f/%.3f spc %.3f/%s mean %.6f range %.2f-%.2f"%(r["shell"],r["rho"],r["sigma_target"],r["sigma_measured"],r["spacing_target"],r["spacing_measured"],r["mean_over_rho"],r["min_over_rho"],r["max_over_rho"])) for r in B["rows"]]
