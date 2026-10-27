# Phase5 (NEW): TT chi2 (ensemble, l=2-2500) and Migkas chi2 vs d for all CMB cases on Migkas-profile backgrounds (summary6.json).
import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
O='/workspace/coepdf/joint_solve/phase5/';R=[r for r in json.load(open(O+'summary6.json')) if 'TT' in r]
fig,ax=plt.subplots(1,2,figsize=(11,4.2))
for grp,mk,lab in (('A_d','o','fixed d, k free (k runs away)'),('K','s','d free, k fixed (phase4b profile)'),('J_','^','fixed d and k')):
    S=[r for r in R if r['tag'].startswith(grp)]
    if not S:continue
    x=[max(r['d'],1e-3) for r in S];ax[0].scatter(x,[r['TT'] for r in S],marker=mk,label=lab);ax[1].scatter(x,[r['mig_chi2'] for r in S],marker=mk,label=lab)
    for r,xx in zip(S,x):ax[0].annotate('k=%.2g'%r['k'],(xx,r['TT']),fontsize=6)
for a,y in zip(ax,('TT $\\chi^2$ ($\\ell$=2-2500, ensemble)','Migkas $\\chi^2$')):a.set_xscale('log');a.set_xlabel('d (Gpc)');a.set_ylabel(y);a.legend(fontsize=7)
plt.suptitle('<H0>=68.64, age 13.787 Gyr fixed; CMB background from the Migkas density profile');plt.tight_layout();plt.savefig(O+'chi2_vs_d6_plot.png',dpi=130)
