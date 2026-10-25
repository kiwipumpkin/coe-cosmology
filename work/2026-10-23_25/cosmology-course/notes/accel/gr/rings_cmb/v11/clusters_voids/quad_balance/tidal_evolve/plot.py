import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=json.load(open('te.json'));fig,ax=plt.subplots(1,3,figsize=(20,6))
for k,c in (('tidal','C3'),('self','C0'),('mix50','C2'),('both_raw','k')):
    s=[0]+[r['rms_disp_Gpc'] for r in d[k]];m2=[d['start']['per_m'][1]]+[r['m2'] for r in d[k]];t=[d['start']['total']]+[r['total'] for r in d[k]];m1=[d['start']['per_m'][0]]+[r['per_m'][0] for r in d[k]]
    lab={'tidal':'CoE tidal only','self':'web self-gravity only','mix50':'50/50 mix','both_raw':'both (CoE 1e54 kg + self)'}[k]
    ax[0].plot(s,m2,'o-',c=c,label=lab);ax[1].plot(s,t,'o-',c=c,label=lab);ax[2].plot(s,m1,'o-',c=c,label=lab)
for a,t in zip(ax,('ring quadrupole m=2 (× CMB)','total ring residual (× CMB)','m=1 (× CMB)')):a.set_xscale('symlog',linthresh=0.01);a.set_xlabel('rms Zel\'dovich displacement (Gpc)');a.set_title(t);a.legend(fontsize=9);a.axvline(0.1,ls=':',c='gray')
plt.savefig('tidal_evolve.png',dpi=110,bbox_inches='tight')
