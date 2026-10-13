"""Schematic: flat density -> seed -> potential dipole -> drainage -> clump+void, on one thin shell (shell_fig.py style).
delta(theta) = seed profile * linear growth a/a_rec; crest seed +1.5e-3, trough seed -2.5e-3; Phi from 1D Poisson Phi''=delta (zero mean)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size":16,"axes.labelsize":16,"xtick.labelsize":14,"ytick.labelsize":14})
SPAN=45; th=np.linspace(0,SPAN,3000); R=3.2; Dl=0.12
G=lambda c,w: np.exp(-0.5*((th-c)/w)**2)
prof_c=G(15,3.5); prof_v=G(31,5.0)
def delta(a): g=a*1100; return 1.5e-3*g*prof_c-2.5e-3*g*prof_v
def phi(d):
    k=np.fft.rfftfreq(len(th),th[1]-th[0])*2*np.pi; F=np.fft.rfft(d-d.mean()); F[1:]/=-(k[1:]**2); F[0]=0; return np.fft.irfft(F,len(th))
ST=[(1/1100,"1. Flat density (recombination)","seed ≈ 10⁻³: invisible, density flat"),
    (0.05,"2. Slight perturbation","crest +, trough − ; shell thickness ∝ density"),
    (0.3,"3. Potential dipole: matter drains","flow runs down the potential gradient, trough → crest"),
    (1.0,"4. Clump and void","crest reaches δc = 1.686, trough δv = −2.81")]
for n,(a,title,sub) in enumerate(ST,1):
    d=delta(a); rho=np.clip(1+d,0.05,None); p=phi(d)
    W,HW,HG=8,5.0,1.9; H=HW+2*HG+2.4; fig=plt.figure(figsize=(W,H))
    fig.text(0.5,1-0.25/H,title,ha="center",fontsize=18,weight="bold"); fig.text(0.5,1-0.62/H,sub+"  (schematic)",ha="center",fontsize=14)
    aw=fig.add_axes([0.04,(2*HG+1.2)/H,0.92,HW/H]); aw.axis("off"); Rm=R*1.08; aw.set_xlim(-0.05,Rm+0.05); aw.set_ylim(-0.25,Rm*np.sin(np.radians(SPAN))+0.05)
    t=np.radians(th); aw.plot([0,Rm],[0,0],color="0.6",lw=1); aw.plot([0,Rm*np.cos(t[-1])],[0,Rm*np.sin(t[-1])],color="0.6",lw=1)
    ro=R+0.5*Dl*rho; ri=R-0.5*Dl*rho
    aw.fill(np.r_[ro*np.cos(t),(ri*np.cos(t))[::-1]],np.r_[ro*np.sin(t),(ri*np.sin(t))[::-1]],color="#08306b",lw=0)
    for c,lab,col in ((15,"crest",  "#b2182b"),(31,"trough","#2166ac")):
        if n>1: aw.text((R+0.32)*np.cos(np.radians(c)),(R+0.32)*np.sin(np.radians(c)),{"crest":"clump" if n==4 else "crest","trough":"void" if n==4 else "trough"}[lab],color=col,fontsize=15,ha="center")
    if n==3:
        for c0,c1 in ((27,19),): 
            for rr in (R-0.25,):
                aw.annotate("",xy=(rr*np.cos(np.radians(c1)),rr*np.sin(np.radians(c1))),xytext=(rr*np.cos(np.radians(c0)),rr*np.sin(np.radians(c0))),arrowprops=dict(arrowstyle="->",lw=3,color="0.3"))
    aw.plot(0,0,"k*",ms=16); aw.text(0.08,-0.15,"CoE",fontsize=14)
    a1=fig.add_axes([0.17,(HG+1.0)/H,0.78,(HG-0.35)/H]); a1.plot(th,rho,color="#08306b",lw=3); a1.axhline(1,color="k",ls="--",lw=1.5)
    a1.set_xlim(0,SPAN); a1.set_ylim(0,3); a1.set_ylabel("density\nρ/ρ̄"); a1.set_xticklabels([])
    a2=fig.add_axes([0.17,0.55/H,0.78,(HG-0.35)/H])
    if n<4:
        sc=max(np.abs(p).max(),1e-12); a2.plot(th,p/np.abs(phi(delta(1.0))).max(),color="#7f3b08",lw=3); a2.axhline(0,color="k",ls="--",lw=1.5)
        a2.set_ylim(-1.1,1.1); a2.set_ylabel("potential\nΦ (rel.)"); a2.set_xlim(0,SPAN); a2.set_xlabel("angle along shell (deg)")
        if n>1: a2.text(15,-0.95 if n==3 else -0.5,"well",ha="center",fontsize=14,color="#b2182b"); a2.text(31,0.75 if n==3 else 0.35,"hill",ha="center",fontsize=14,color="#2166ac")
    else:
        aa=np.logspace(np.log10(1/1100),0,200)
        a2.semilogx(aa,1.5e-3*1100*aa,color="#b2182b",lw=3); a2.semilogx(aa,-2.5e-3*1100*aa,color="#2166ac",lw=3)
        a2.axhline(1.686,color="#b2182b",ls="--",lw=1.5); a2.axhline(-2.81,color="#2166ac",ls="--",lw=1.5); a2.axhline(0,color="k",ls="--",lw=1)
        a2.set_ylim(-3.4,2.6); a2.set_ylabel("linear δ"); a2.set_xlabel("scale factor (recombination → today)")
        a2.text(1.2e-3,1.9,"δc = 1.686",color="#b2182b",fontsize=14); a2.text(1.2e-3,-2.5,"δv = −2.81",color="#2166ac",fontsize=14)
    fig.savefig(f"flat_to_clump_{n}.png",dpi=150); plt.close(fig)
