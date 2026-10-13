"""Chris's sheet definition (spreadsheet_model/columns.tex, col AF):
A/A_all(th) = | kR n + dt [g(r_P) - g(D x)] | / (kR), P at distance R along n. Point mass at CoE,
strength eps = G M dt /(D^2 k R) = 0.247 (sheet value). D = 58.2, R = 14.3 Gpc. th from the CoE direction."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D,R,eps=58.2,14.3,0.247
ang=np.linspace(0,360,1441); t=np.radians(ang)
n=np.stack([-np.cos(t),np.sin(t)],1)              # ang=0 points at the CoE (CoE at origin, observer at (D,0))
xo=np.array([D,0.]); xp=xo+R*n; rp=np.linalg.norm(xp,axis=1)
g=lambda x,r:-x/r[:,None]**3
dg=(g(xp,rp)-g(xo[None],np.array([D])))*D**2       # in units of GM/D^2
AF_mass=np.linalg.norm(n+eps*dg,axis=1); AF_v=np.ones_like(ang)
mean=np.trapezoid(AF_mass,ang)/360
fig,ax=plt.subplots(figsize=(8.5,4.6))
ax.plot(ang,AF_v,color="0.4",lw=2,ls="--",label="velocity only, V = k r (no mass)")
ax.plot(ang,AF_mass,color="#c0392b",lw=2.2,label=f"with point mass at CoE (GM·dt/D²kR = {eps})")
for a in (90,270): ax.axvline(a,color="0.7",lw=0.8,ls=":")
ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_ylim(0,max(1.1,AF_mass.max()*1.08))
ax.set_xlabel("angle from CoE direction (deg)",fontsize=11); ax.set_ylabel("A / A_all  (speed with gravity / without)",fontsize=11)
ax.set_title("Chris's sheet quantity (col AF) vs angle, D = 58.2, R = 14.3 Gpc",fontsize=12); ax.legend(frameon=False,fontsize=10,loc="lower center")
fig.savefig("migkas_ratio.png",dpi=130,bbox_inches="tight")
i=np.argsort(AF_mass); print("min",AF_mass.min(),"at",sorted(set(np.round(ang[i[:2]]))),"max",AF_mass.max(),"at",ang[np.argmax(AF_mass)],"mean",mean,"AF at 0,180",AF_mass[0],AF_mass[720])
from scipy.signal import argrelmin; print("local minima at",ang[argrelmin(AF_mass)[0]])
