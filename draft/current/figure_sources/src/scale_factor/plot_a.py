import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.integrate import quad
Gyr=977.8  # 1/(km/s/Mpc) in Gyr
def t_of_a(a,H0,Om,OL):
    Ok=1-Om-OL
    f=lambda x: 1/(x*np.sqrt(Om/x**3+Ok/x**2+OL))
    return quad(f,1e-8,a)[0]*Gyr/H0
A=np.linspace(1e-3,3,600)
cases=[("ΛCDM (Planck: H0=67.4, Ωm=0.315, ΩΛ=0.685)",67.4,0.315,0.685,'k','-'),
       ("CoE shell, no mass (Milne, H0=70)",70,0.0,0,'C0','--'),
       ("CoE shell, mass parameter 0.1 (H0=70)",70,0.1,0,'C1','--'),
       ("CoE shell, mass parameter 0.3 (H0=70)",70,0.3,0,'C3','--')]
fig,ax=plt.subplots(figsize=(8,5.5))
rows=[]
for lab,H0,Om,OL,c,ls in cases:
    t=np.array([t_of_a(a,H0,Om,OL) for a in A]); t0=t_of_a(1,H0,Om,OL)
    ax.plot(t,A,color=c,ls=ls,lw=2,label=f"{lab}; age {t0:.1f} Gyr")
    ax.plot(t0,1,'o',color=c)
    rows.append(np.column_stack([t,A]))
ax.axhline(1,color='grey',lw=.5); ax.set_xlim(0,40); ax.set_ylim(0,3)
ax.set_xlabel("time since bang (Gyr)"); ax.set_ylabel("scale factor a (a = 1 today)")
ax.set_title("Scale factor vs time: ΛCDM vs CoE shell (exact GR, E > 0)")
ax.legend(fontsize=8,loc='upper left'); fig.tight_layout(); fig.savefig("scale_factor_vs_time.png",dpi=150)
np.savez("scale_factor_vs_time.npz",a=A,**{f"t_case{i}":r[:,0] for i,r in enumerate(rows)})
for lab,H0,Om,OL,*_ in cases: print(lab, round(t_of_a(1,H0,Om,OL),2))
