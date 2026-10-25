import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
exec(open('cum2.py').read().split("rA=Araw")[0])
d=json.load(open('../l50/curves.json'));g0=np.array(d['grid']);A0=np.array(d['A0']);B0=np.array(d['B0']);AFr=np.array(d['AF']);BLr=np.array(d['BL'])
# sample reconstruction on the 71 bin centres
A0l=periodic_interp(g0,A0,l);B0l=periodic_interp(g0,B0,l);AFl=periodic_interp(g0,AFr,l);BLl=periodic_interp(g0,BLr,l)
rA=A0l-AFl;rB=B0l-BLl;base=(rA.copy(),rB.copy());steps2=[]
for n in range(1,5):
  h=np.column_stack([np.cos(n*th),np.sin(n*th)]);X=np.vstack([np.column_stack([o,Z,h]),np.column_stack([Z,o,h])]);s=np.r_[sA,sB]
  p=np.linalg.lstsq(X/s[:,None],np.r_[rA,rB]/s,rcond=None)[0];steps2.append((n,p));rA=rA-X[:71]@p;rB=rB-X[71:]@p
for n,p in steps2:
  C=np.hypot(p[2],p[3]);ph=np.degrees(np.arctan2(p[3],p[2]))%360;print(n,f'C={C:.3f}',sorted(round((ph/n+360*j/n)%360,1) for j in range(n)))
run=lambda r:np.array([np.sqrt(np.mean(r[np.abs((l-c+180)%360-180)<=15]**2)) for c in l])
plt.rcParams.update({'font.size':15});fig,axs=plt.subplots(2,1,figsize=(14,10),facecolor='white',sharex=True)
for a,nm,b0,fn,c,sg in ((axs[0],'A',base[0],rA,'tab:blue',sA),(axs[1],'B',base[1],rB,'tab:red',sB)):
  a.axvspan(15,80,color='gold',alpha=.18);a.axvspan(260,335,color='lightsteelblue',alpha=.35);a.axhline(0,color='grey',lw=1)
  a.plot(l,b0,'o-',color='k',lw=2,ms=4,label=f'data − CoE model (no harmonics), RMS {np.sqrt(np.mean(b0**2)):.3f}')
  a.plot(l,fn,'s--',color=c,lw=2,ms=4,label=f'after shared n=1–4, RMS {np.sqrt(np.mean(fn**2)):.3f}')
  a.plot(l,run(b0),color='k',lw=1.2,ls=':',label='running RMS ±15° (no harmonics)');a.plot(l,run(fn),color=c,lw=1.2,ls=':',label='running RMS ±15° (after n=1–4)')
  a.set_ylabel(f'{nm} residual');a.set_ylim(-0.15,0.15);a.legend(loc='upper center',bbox_to_anchor=(.5,-.04),ncol=2,frameon=False,fontsize=11.5)
axs[1].set_xlim(0,360);axs[1].set_xticks(range(0,361,30));axs[1].set_xlabel('Galactic longitude ℓ [deg]')
fig.suptitle('Residuals vs angle — reconstruction (313 clusters) vs CoE model (A: AF_smooth, B: sheet BL)')
fig.tight_layout(h_pad=3);fig.savefig('recon_resid_vs_angle.png',dpi=120,facecolor='white')
