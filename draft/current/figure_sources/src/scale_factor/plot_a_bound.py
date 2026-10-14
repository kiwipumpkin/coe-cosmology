import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
Gyr=977.8
def run(H0,Om,OL):
    h=H0/Gyr  # 1/Gyr
    f=lambda t,y:[y[1], -0.5*h**2*Om/y[0]**2 + h**2*OL*y[0]]
    ev=lambda t,y:y[0]-1e-3; ev.terminal=True
    b=solve_ivp(f,[0,-60],[1,h],events=ev,max_step=0.01,rtol=1e-9)
    t0=-b.t[-1]
    fw=solve_ivp(f,[0,400],[1,h],events=ev,max_step=0.01,rtol=1e-9)
    t=np.concatenate([b.t[::-1],fw.t[1:]])+t0; a=np.concatenate([b.y[0][::-1],fw.y[0][1:]])
    return t,a,t0
cases=[("ΛCDM flat (Planck Ωm=0.315, ΩΛ=0.685)",67.4,0.315,0.685,'k','-'),
 ("ΛCDM closed, E<0 (Ωm=0.315, ΩΛ=0.735, Ωk=−0.05)",67.4,0.315,0.735,'grey','-'),
 ("ΛCDM closed, E<0, heavy (Ωm=2.0, ΩΛ=0.685)",67.4,2.0,0.685,'m','-'),
 ("CoE shell, E>0 (mass 0.3)",70,0.3,0,'C3','--'),
 ("CoE shell, E=0 critical (mass 1)",70,1.0,0,'C1','--'),
 ("CoE shell, E<0 bound (mass 1.5)",70,1.5,0,'C2','--'),
 ("CoE shell, E<0 bound (mass 3)",70,3.0,0,'C4','--')]
fig,ax=plt.subplots(figsize=(9,6)); out={}
for i,(lab,H0,Om,OL,c,ls) in enumerate(cases):
    t,a,t0=run(H0,Om,OL); k=np.argmax(a)
    note=f"; turns around at {t[k]:.0f} Gyr, crunch {t[-1]:.0f} Gyr" if a[-1]<0.01 else ""
    ax.plot(t,a,color=c,ls=ls,lw=2,label=f"{lab}; age {t0:.1f} Gyr{note}"); ax.plot(t0,1,'o',color=c)
    out[f"t{i}"]=t; out[f"a{i}"]=a; print(lab,round(t0,2),note)
ax.axhline(1,color='grey',lw=.5); ax.set_xlim(0,160); ax.set_ylim(0,3.5)
ax.set_xlabel("time since bang (Gyr)"); ax.set_ylabel("scale factor a (a = 1 today)")
ax.set_title("Scale factor vs time: E > 0, E = 0 and E < 0 cases")
ax.legend(fontsize=7.5,loc='lower right'); fig.tight_layout(); fig.savefig("scale_factor_bound.png",dpi=150)
np.savez("scale_factor_bound.npz",**out)
