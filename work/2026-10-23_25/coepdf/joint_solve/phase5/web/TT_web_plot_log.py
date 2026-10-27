# Phase5 web test: TT of Route A eps=0.01 + web compensation vs Planck, l=2..2500. LOG=1 -> log l axis (TT_web_plot_log.py is the same code).
import json,os,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
O='/workspace/coepdf/joint_solve/phase5/';b=json.load(open(O+'web/tt_web.json'));t=b['tag'];z=np.load(O+f'p6_{t}.npz');LOG=os.environ.get('LOG','0')=='1'
D,Drc,Lt,Pf=z['D'],z['Drc'],z['Lt'],z['Pf'];l=np.arange(len(D));m=(l>=2)&(l<=2500)
fig,ax=plt.subplots(figsize=(10,5));ax.errorbar(Pf[:,0],Pf[:,1],yerr=[Pf[:,2],Pf[:,3]],fmt='k.',ms=2,elinewidth=0.3,alpha=0.5,label='Planck 2018 TT data');ax.plot(np.arange(len(Lt))[2:2501],Lt[2:2501],'g-',lw=0.8,label='Planck best-fit LCDM')
ax.plot(l[m],D[m],'r-',lw=1.2,label='model ensemble mean');ax.plot(l[m],Drc[m],'b-',lw=0.6,alpha=0.6,label='model seed-42 sky (compensated)')
ax.set_xlim(2,2500);ax.set_xlabel('$\\ell$');ax.set_ylabel('$D_\\ell$ ($\\mu K^2$)')
if LOG:ax.set_xscale('log')
ax.set_title(f"Route A eps=0.01 + cluster/void web: TT chi2={b['TT']:.0f} (l2-30 {b['TT_2_30']:.0f})",fontsize=9)
ax.legend();plt.tight_layout();plt.savefig(O+('web/TT_web_plot_log.png' if LOG else 'web/TT_web_plot.png'),dpi=130)
