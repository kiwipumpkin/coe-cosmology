import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import brentq
c=1.0; t=3; d=3.0+t
def bands(K):
    L=np.linspace(d+c*t,d-c*t,2*K*t+1); ro,ri=L[:-1],L[1:]
    return ro,ri,(ro+ri)/2,(4/3)*np.pi*(ro**3-ri**3)
f1,a1=plt.subplots(figsize=(10,6.5)); f2,a2=plt.subplots(figsize=(10,6.5))
for K,ls,lw in ((1,"-",2.4),(2,"--",1.6)):
    ro,ri,m,V=bands(K); x=np.r_[ri[::-1],ro[0]]
    p=brentq(lambda p:(m[-1]**p/V[-1])/(m[0]**p/V[0])-1.03,0.5,4)
    for q,lab,col in ((np.ones_like(m),"equal mass per shell","tab:red"),(m**p,f"3% lock (mass ∝ r^{p:.3f})","tab:purple")):
        qn=q/q[-1]; rho=q/V; rn=rho/rho[-1]
        a1.step(x,np.r_[qn[::-1],qn[0]],where="post",color=col,ls=ls,lw=lw,label=f"{lab}, {6*K} shells")
        a2.step(x,np.r_[rn[::-1],rn[0]],where="post",color=col,ls=ls,lw=lw,label=f"{lab}, {6*K} shells")
        print(K,lab,"mass out/in",qn[0],"density out/in",rn[0])
    assert abs(1/((m[0]**p/V[0])/(m[-1]**p/V[-1]))-1.03)<1e-10
for a,yl,tt in ((a1,"mass+light per shell ÷ innermost shell","Mass+light per shell, t = 3"),
                (a2,"density ÷ innermost band (log scale)","Density, t = 3")):
    a.axvline(3,color="k",ls=":"); a.axvline(9,color="k",ls=":")
    a.set_xlabel("distance from CoE (c·Δt)"); a.set_ylabel(yl); a.set_title(tt+"   (solid = 6 shells, dashed = 12 shells)")
    a.grid(alpha=.3,which="both"); a.set_xlim(2.8,9.2)
a2.set_yscale("log")
for a in (a1,a2):
    y=a.get_ylim()[1]; a.text(3,y," D−R",va="top"); a.text(9,y,"D+R ",va="top",ha="right")
a1.legend(loc="upper left",bbox_to_anchor=(0.05,0.9)); a2.legend(loc="lower left")
f1.tight_layout(); f1.savefig("mass_per_shell_both.png",dpi=180); f2.tight_layout(); f2.savefig("density_both.png",dpi=180)
