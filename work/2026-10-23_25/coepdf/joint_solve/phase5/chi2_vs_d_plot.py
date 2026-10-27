# Phase5 (NEW): TT chi2 (ensemble mean; l=2-2500 and l=2-30) and Migkas chi2 vs d, from summary5.json.
import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
O='/workspace/coepdf/joint_solve/phase5/';R=[r for r in json.load(open(O+'summary5.json')) if 'TT_mean' in r];d=np.array([max(r['d'],1e-3) for r in R])
fig,ax=plt.subplots(3,1,figsize=(8,9),sharex=True)
ax[0].plot(d,[r['TT_mean'] for r in R],'ro-',label='ensemble mean');ax[0].plot(d,[r['TT_comp'] for r in R],'s--',color='orange',label='seed-42 compensated');ax[0].set_ylabel(r'TT $\chi^2$ ($\ell$=2-2500)');ax[0].legend()
ax[1].plot(d,[r['TT_mean_2_30'] for r in R],'ro-');ax[1].plot(d,[r['TT_comp_2_30'] for r in R],'s--',color='orange');ax[1].set_ylabel(r'TT $\chi^2$ ($\ell$=2-30)')
ax[2].plot(d,[r['mig_chi2'] for r in R],'ko-');ax[2].set_ylabel(r'Migkas $\chi^2$ (n=71)');ax[2].set_xlabel('d [Gpc]');ax[2].set_xscale('symlog',linthresh=0.5)
for a in ax:a.grid(alpha=0.3)
ax[0].set_title(r'Phase 5: d scan, $\langle H_0\rangle$=68.64, age 13.787 Gyr')
plt.tight_layout();plt.savefig(O+'chi2_vs_d_plot.png',dpi=140)
