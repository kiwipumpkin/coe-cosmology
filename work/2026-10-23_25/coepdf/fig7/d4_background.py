# REGENERATED for rev7 (replaces rev5/6 radiation+matter d4_radmat with t0=13.8 Gyr): the background actually used by the solve.
# LTB famP observer-shell history from los_flow/ltb.npz; consistency R = eta_rec * D from los_flow_R/tab.npz; H0, age from p3_background_relations.json.
import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
z=np.load('/workspace/redweb2/realweb/los_flow/ltb.npz');tb=np.load('/workspace/redweb2/realweb/los_flow_R/tab.npz')
B=json.load(open('/workspace/coepdf/joint_solve/phase2/p3_background_relations.json'));b=B['by_D']['34.96'];D=34.96;Gyr=D*3.26156
t=z['T']*Gyr;fig,ax=plt.subplots(1,3,figsize=(13,4.2))
ax[0].loglog(t,z['a_t'],'k',lw=2,label='a(t), observer shell');ax[0].loglog(t,(t/t[-1])**(2/3),'g:',label=r'$t^{2/3}$');ax[0].set_xlabel('t [Gyr]');ax[0].set_ylabel('a')
ax[1].loglog(t,z['H']/z['H'][-1],'b',lw=2,label=r'$H/H_0$');ax[1].loglog(t,z['Omega'],'r',lw=1.5,label=r'$\Omega$ (local)');ax[1].set_xlabel('t [Gyr]')
ax[1].set_title(f"$H_0$={b['H0_ltb_kms_Mpc']:.1f} km/s/Mpc, age={b['age_ltb_Gyr']:.3f} Gyr",fontsize=10)
Ds=[float(k) for k in B['by_D']];ax[2].plot(Ds,[B['by_D'][k]['R_ls_ltb_Gpc'] for k in B['by_D']],'ko-',label=r'$R_{ls}=\eta_{\rm rec}D$')
ax[2].axvline(D,color='m',ls='--',label=f'chosen D={D} Gpc, R={b["R_ls_ltb_Gpc"]:.4f} Gpc');ax[2].set_xlabel('D [Gpc]');ax[2].set_ylabel('R [Gpc]')
for a in ax:a.legend(fontsize=8);a.grid(alpha=.3)
fig.tight_layout();fig.savefig('d4_background.png',dpi=140,facecolor='white')
