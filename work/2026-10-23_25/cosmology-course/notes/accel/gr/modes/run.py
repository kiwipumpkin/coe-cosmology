import json,numpy as np;from modes import *
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':15})
TH=np.arange(360.);names={1:'observer & sphere at rest',2:'observer moving, sphere at rest',3:'observer + sphere moving together',4:'together, accelerating along motion',5:'together, accelerating against motion'}
res=[]
for mode in range(1,6):
  for ex in (False,True):
    fig,ax=plt.subplots(figsize=(12,7));tag='expanding' if ex else 'static'
    for v,a,c in [(0.01,0.01,'C0'),(0.1,0.01,'C3'),(0.1,0.02,'C2')]:
        if mode in(1,) and v!=0.01:continue
        if mode in(2,3) and a!=0.01:continue
        s,d=stats(ring(**build(mode,v,a,ex)));res.append(dict(mode=mode,sphere=tag,v=v,a=(a if mode>3 else 0),**s))
        lab=('' if mode==1 else f'v={v}c'+(f', a={a} c²/R' if mode>3 else ''))+f": RMS={s['rms']:.2e}, range={s['range']:.2e}, dipole={s['dipole']:.2e}, quad={s['quad']:.2e}"
        ax.plot(TH,d,c,lw=2.5,label=lab)
    ax.set_xlabel('sight angle from direction of motion, observer frame (deg)');ax.set_ylabel('ΔT/T');ax.set_title(f'Mode {mode}: {names[mode]} ({tag} sphere)',fontsize=13)
    ax.legend(loc='lower center',bbox_to_anchor=(.5,1.06),frameon=False,fontsize=11);fig.tight_layout();fig.savefig(f'mode{mode}_{tag}.png',dpi=105);plt.close()
json.dump(res,open('modes.json','w'),indent=1)
for r in res:print(r['mode'],r['sphere'],r['v'],r['a'],'rms %.3e range %.3e dip %.3e(%+d) quad %.3e mono %.5f'%(r['rms'],r['range'],r['dipole'],r['dipole_sign'],r['quad'],r['monopole_T_over_Temit']))
