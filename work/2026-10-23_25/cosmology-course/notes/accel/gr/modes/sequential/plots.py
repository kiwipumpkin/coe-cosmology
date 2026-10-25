import pickle,json,numpy as np;from seq import *
from scipy.optimize import brentq
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':13})
R=pickle.load(open('seq_main.pkl','rb'));R.update(pickle.load(open('seq_p092.pkl','rb')));TH=np.arange(0,360,2.0)
for p in (1.0,1.08,0.92):
    fig,ax=plt.subplots(3,1,figsize=(12,13),sharex=True)
    for (pp,(rule,ar)),out in R.items():
        if pp!=p:continue
        T=np.array([1/o['opz'] if o else np.nan for o in out]);d=T/np.nanmean(T)-1
        lab=f"recomb. scale {ar}: RMS={np.sqrt(np.nanmean(d**2)):.2e}, range={np.nanmax(d)-np.nanmin(d):.2e}"
        ax[0].plot(TH,d,lw=2,label=lab);ax[1].plot(TH,[o['te'] if o else np.nan for o in out],lw=2);ax[2].plot(TH,[o['travel'] if o else np.nan for o in out],lw=2)
    ax[0].set_ylabel('ΔT/T');ax[1].set_ylabel('emission time (D/c)');ax[2].set_ylabel('light travel time = path length (D/c)');ax[2].set_xlabel('sight angle (deg; 0 = away from CoE, 180 = toward CoE)')
    ax[0].legend(loc='lower center',bbox_to_anchor=(.5,1.02),frameon=False,fontsize=11);ax[0].set_title(f'velocity law u ∝ r^{p}'+(' (Milne)' if p==1 else ''),pad=95)
    fig.tight_layout();fig.savefig(f'seq_p{p}.png',dpi=100);plt.close()
# 3000 K surface: areal radius where local scale = 0.01, vs time
fig,ax=plt.subplots(figsize=(11,7));out={}
for p in (1.08,0.92):
    u=U(p);ts=np.linspace(0.02,0.9,40);Rs=[]
    for t in ts:
        rg=np.geomspace(1e-3,200,4000);a=np.array([np.cbrt(u.f(r,t)[0]**2*u.f(r,t)[1]/r**2) if t>u.tB(r) else np.nan for r in rg])-0.01
        s=np.where(np.isfinite(a[:-1])&np.isfinite(a[1:])&(np.sign(a[:-1])!=np.sign(a[1:])))[0]
        Rs.append(float(u.f(rg[s[0]],t)[0]) if len(s) else np.nan)
    out[p]=dict(t=ts.tolist(),R3000=Rs);ax.plot(ts,Rs,lw=2.5,label=f'u ∝ r^{p}')
ax.set_xlabel('time (D/c, today = %.3f)'%T1);ax.set_ylabel('areal radius of 3000 K surface (D), recomb. scale 0.01');ax.set_yscale('log');ax.legend(loc='lower center',bbox_to_anchor=(.5,1),frameon=False,ncol=2)
fig.tight_layout();fig.savefig('R3000_vs_t.png',dpi=100);json.dump(out,open('R3000.json','w'))
for p in out: print(p,[ (round(t,2),None if not np.isfinite(r) else round(r,3)) for t,r in zip(out[p]['t'][::6],out[p]['R3000'][::6])])
