import json, numpy as np, tcore as T, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
R=json.load(open("solve.json")); Rg=14.3; col={"0.03-0.07":"C1","z>0.07":"C3","all":"k"}; nm={"0.03-0.07":"orange 0.03-0.07","z>0.07":"red z>0.07","all":"black all"}
fig,ax=plt.subplots(figsize=(11,6.8)); Dg=np.geomspace(14.75,200,400); dmax={}
for s in col:
    rr=sorted([r for r in R if r["slice"]==s and r["mode"]=="median" and not r["cone"] and r["l0_branch"]==0],key=lambda r:r["D"])
    D=np.array([r["D"] for r in rr]); f=np.array([r["mcc"] for r in rr]); rc=np.array([T.rhoc(r["H_E"]) for r in rr])
    fg=np.interp(Dg,D,f); rcg=np.interp(Dg,D,rc)
    Mc=fg*4/3*np.pi*rcg*(Dg*T.GPCm)**3; Ms=0.049*rcg*8*np.pi*Rg*(Dg**2+Rg**2/3)*T.GPCm**3
    ax.plot(Dg,Mc,'-',color=col[s],lw=2.2,label=nm[s]+": core mass"); ax.plot(Dg,Ms,'--',color=col[s],lw=1.4,label=nm[s]+": shell mass")
    ratio=fg*(Dg/(Dg-Rg))**3; dmax[s]=Dg[np.where(ratio>=0.049)[0].max()]
ax.axvline(14.7,color='C2',lw=2); ax.text(15.1,3.5e56,'D min\n14.7 Gpc\n(see out\nto R)',color='C2',fontsize=9.5,va='top')
for s in col: ax.axvline(dmax[s],color=col[s],lw=1.6,ls=':')
ax.axvspan(14.7,dmax["all"],color='C2',alpha=.08); ax.axvspan(dmax["all"],220,color='C0',alpha=.06)
ax.text(15.1,1.6e51,'red\n%.1f'%dmax["z>0.07"],color='C3',fontsize=9.5)
ax.text(dmax["all"]*1.03,2.2e51,'black %.1f'%dmax["all"],color='k',fontsize=9.5)
ax.text(dmax["0.03-0.07"]*1.03,6e51,'orange %.1f'%dmax["0.03-0.07"],color='C1',fontsize=9.5)
ax.text(15.1,1.2e53,'core\ndenser\nthan shell\n(D max)',color='C2',fontsize=9.5)
ax.text(26,3.5e56,'core less dense than shell (Sedov-like thin interior):\nthe dotted lines become D min, no D max',color='C0',fontsize=10,va='top')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(12,220); ax.set_ylim(8e50,6e56)
ax.set_xticks([14.7,20,30,49,94,200]); ax.set_xticklabels(['14.7','20','30','49','94','200']); ax.xaxis.set_minor_formatter(NullFormatter())
ax.set_xlabel('CoE distance D (Gpc)'); ax.set_ylabel('mass (kg)'); ax.grid(alpha=.3,which='both'); ax.legend(fontsize=9,loc='lower right')
ax.set_title('Core and shell mass vs D. Dotted lines: core density (spread to D - R) equals the shell, 0.049 of critical')
plt.tight_layout(); plt.savefig('mass_vs_D_limits_v6.png',dpi=120)
