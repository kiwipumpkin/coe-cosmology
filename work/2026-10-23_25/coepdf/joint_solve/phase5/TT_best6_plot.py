# Phase5 (NEW): TT spectrum of the best joint (Migkas+TT) case vs Planck, linear l axis. Set LOG=1 for log l axis (TT_best6_plot_log.py).
import json,os,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
O='/workspace/coepdf/joint_solve/phase5/';b=json.load(open(O+'best6.json'));t=b['tag'];z=np.load(O+f'p6_{t}.npz');LOG=os.environ.get('LOG','0')=='1'
D,Drc,Lt,Pf=z['D'],z['Drc'],z['Lt'],z['Pf'];l=np.arange(len(D));m=(l>=2)&(l<=2500)
fig,ax=plt.subplots(figsize=(10,5));ax.plot(Lt,Pf,'k.',ms=2,label='Planck 2018 TT')
ax.plot(l[m],D[m],'r-',lw=1.2,label='model ensemble mean');ax.plot(l[m],Drc[m],'b-',lw=0.6,alpha=0.6,label='model seed-42 sky (compensated)')
ax.set_xlim(2,2500);ax.set_xlabel('$\\ell$');ax.set_ylabel('$D_\\ell$ ($\\mu K^2$)')
if LOG:ax.set_xscale('log')
ax.set_title(f"{t}: d={b['d']:.3g} Gpc, k={b['k']:.3g}, eta={b['eta']:.3g}, D={b['D']:.2f} Gpc; TT chi2={b['TT']:.0f} (l2-30 {b['TT_2_30']:.0f}), Migkas chi2={b['mig_chi2']:.3f}",fontsize=8)
ax.legend();plt.tight_layout();plt.savefig(O+('TT_best6_plot_log.png' if LOG else 'TT_best6_plot.png'),dpi=130)
