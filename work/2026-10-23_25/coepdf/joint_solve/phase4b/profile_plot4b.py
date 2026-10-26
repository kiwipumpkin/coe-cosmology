# Phase4b (NEW): plot profile chi2 vs k from profile_k4b.json (A: d free, B: d=14.30 Gpc).
import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
O='/workspace/coepdf/joint_solve/phase4b/';R=json.load(open(O+'profile_k4b.json'))
fig,ax=plt.subplots(figsize=(6,4))
for m,lab in (('A','d free (down to 1e-3 Gpc)'),('B','d = 14.30 Gpc')):
    r=sorted([(lk,f) for mm,lk,f,_ in R if mm==m])
    if r:x,y=np.array(r).T;ax.plot(10**x,y,'o-',label=lab)
ax.axvline(1,ls=':',c='k');ax.set_xscale('log');ax.set_xlabel('k = gradient / CatWISE');ax.set_ylabel(r'Migkas $\chi^2$ + age penalty (n=71)')
ax.set_title(r'Profile in k, $\langle H_0\rangle$=68.64, age 13.787 Gyr (phase4b)');ax.legend();fig.tight_layout();fig.savefig(O+'profile_plot4b.png',dpi=150)
