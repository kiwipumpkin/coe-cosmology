import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, json
from scipy.optimize import brentq
c=2.998e8; G=6.674e-11; GPC=3.0857e25; MPC=GPC/1e3; GYR=3.156e16; KM=1e3
plt.rcParams.update({"font.size":17,"axes.labelsize":18,"axes.titlesize":18,"legend.fontsize":15,"xtick.labelsize":16,"ytick.labelsize":16})
k=160*KM/MPC; T=7.86*GYR; D=3.2*GPC; M1=6.08e51
def run(r0,M,N=4000):
    GM=G*M; h=T/N; r=r0; v=k*r0; ts=[0];rs=[r];vs=[v]
    f=lambda r:-GM/r**2
    for i in range(N):
        a1,b1=v,f(r);a2,b2=v+.5*h*b1,f(r+.5*h*a1);a3,b3=v+.5*h*b2,f(r+.5*h*a2);a4,b4=v+h*b3,f(r+h*a3)
        r+=h/6*(a1+2*a2+2*a3+a4); v+=h/6*(b1+2*b2+2*b3+b4); ts.append((i+1)*h);rs.append(r);vs.append(v)
    return np.array(ts),np.array(rs),np.array(vs)
out={}; cols={1:"#1b9e77",4:"#d95f02",10:"#7570b3"}
lam=np.logspace(0,np.log10(3000),300)*MPC
fig1,a1=plt.subplots(figsize=(8,6)); fig2,a2=plt.subplots(figsize=(8,6))
for m in (1,4,10):
    M=m*M1; r0=brentq(lambda x:run(x,M,800)[1][-1]-D,0.05*D,D)
    t,r,v=run(r0,M); g=G*M/r**2; acc=-np.gradient(v,t)
    a1.plot(t/GYR,g*1e11,color=cols[m],lw=3.5,label="%d×"%m); a1.plot(t/GYR,acc*1e11,"k:",lw=1.5)
    # fixed angular mode: physical k(t)= 2pi/lam * D/r(t)
    Nf=np.array([np.trapezoid(np.sqrt(2*np.pi/L*D/r*g),t) for L in lam])
    a2.plot(lam/MPC,Nf,color=cols[m],lw=3.5,label="%d×"%m)
    sig_now=np.sqrt(2*np.pi/(100*MPC)*g[-1])
    out[m]=dict(r0_Gpc=r0/GPC,v_now_kms=v[-1]/KM,g_start=g[0],g_now=g[-1],decel_check_maxrelerr=float(np.max(np.abs(acc[5:-5]/g[5:-5]-1))),
      efolds={str(L):float(np.interp(L,lam/MPC,Nf)) for L in (10,30,50,100,150,300,1000)},
      lam_1efold_Mpc=float(np.interp(1,Nf[::-1],lam[::-1]/MPC)),lam_3efold_Mpc=float(np.interp(3,Nf[::-1],lam[::-1]/MPC)),
      sigma100_now_perGyr=sig_now*GYR)
a1.set_yscale("log");a1.set_xlabel("time since start (Gyr)");a1.set_ylabel("g_eff = GM/r²  (10⁻¹¹ m/s²)")
a1.set_title("Deceleration of our shell (dots: −d²r/dt²)");a1.legend(frameon=False);fig1.tight_layout();fig1.savefig("p1_geff.png",dpi=110)
a2.set_xscale("log");a2.axhline(1,color="k",ls="--",lw=2);a2.text(1.3,1.1,"1 e-fold",fontsize=15)
for L,nm in ((40,"voids"),(100,"walls"),(150,"BAO")): a2.axvline(L,color="gray",lw=1.5,ls=":");a2.text(L*1.05,{"voids":100,"walls":90,"BAO":80}[nm],nm,fontsize=14,color="gray")
a2.set_xlabel("wavelength today (Mpc)");a2.set_ylabel("RT e-folds  ∫σ dt");a2.set_title("Rayleigh–Taylor growth, A = 1");a2.legend(frameon=False,loc="center right");fig2.tight_layout();fig2.savefig("p2_efolds.png",dpi=110)
# p3 Vishniac spacing
dR=np.linspace(0.01,0.2,100); f,ax=plt.subplots(figsize=(8,6))
for R,nm,ls in ((D,"shell at D = 3.2 Gpc","-"),(14.3*GPC,"shell R = 14.3 Gpc","--")):
    ax.fill_between(dR,2*dR*R/MPC,5*dR*R/MPC,alpha=.25); ax.plot(dR,3*dR*R/MPC,ls,lw=3.5,label=nm+" (λ≈3Δ)")
for L,nm in ((40,"voids 30–50"),(100,"walls ~100"),(150,"BAO 150")): ax.axhline(L,color="gray",ls=":",lw=2); ax.text(0.03 if L==40 else 0.15,L*1.08,nm,fontsize=14,color="gray")
ax.set_yscale("log");ax.set_xlabel("shell thickness Δ/R");ax.set_ylabel("Vishniac spacing λ (Mpc)");ax.set_title("Predicted spacing (band: λ = 2Δ–5Δ)");ax.legend(frameon=False,loc="upper left",fontsize=14);f.tight_layout();f.savefig("p3_spacing.png",dpi=110)
for R in (3.2,14.3): out["lmode_%g"%R]={str(x):[3*x*R*1e3, 2*np.pi/(3*x)] for x in (0.01,0.03,0.1,0.2)}
json.dump(out,open("numbers.json","w"),indent=1); print(json.dumps(out,indent=1))
