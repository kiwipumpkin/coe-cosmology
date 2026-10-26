# REGENERATED for rev7: volume vs time in the consistent LTB famP background (D=34.96 Gpc, R=R_ls=0.18957 Gpc).
# a(t) and H(t) are the observer-shell values from los_flow/ltb.npz (dimensionless, time in D/c); converted with D/c = D[Gpc]*3.2616 Gyr.
import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt,json
z=np.load('/workspace/redweb2/realweb/los_flow/ltb.npz');B=json.load(open('/workspace/coepdf/joint_solve/phase2/p3_background_relations.json'))
D=34.96;R=B['by_D']['34.96']['R_ls_ltb_Gpc'];Gyr=D*3.26156;t=z['T']*Gyr;a=z['a_t'];trec=0.0054225  # eta_rec*D = R
V=lambda r:4/3*np.pi*(r*a)**3
fig,ax=plt.subplots(3,1,figsize=(8,9.5),sharex=True)
ax[0].plot(t,V(D),color='tab:blue',lw=2,label=f'D sphere (D={D} Gpc today)');ax[0].plot(t,V(D+R),color='tab:green',lw=1,label='D+R sphere');ax[0].plot(t,V(R),'tab:red',label=f'R sphere (R=R_ls={R:.3f} Gpc today)')
ax[0].set_ylabel('volume [Gpc³]');ax[0].set_title('Volume vs time on the observer shell, LTB famP background')
ax[1].plot(t,100*V(R)/V(D+R),color='tab:red',label='V(R)/V(D+R)');ax[1].set_ylabel('percent [%]');r=100*V(R)/V(D+R);ax[1].set_ylim(r.min()*0.9,r.max()*1.1);ax[1].ticklabel_format(axis='y',style='sci',scilimits=(0,0),useOffset=False)
H=z['H']/Gyr;ax[2].plot(t,3*H,color='tab:blue',label=r'$\dot V/V=3H$ (ltb.npz)');ax[2].plot(t,2/t,':',color='.4',label='2/t (matter)');ax[2].set_ylabel(r'$\dot V/V$ [Gyr$^{-1}$]');ax[2].set_xlabel('t [Gyr]')
for x in ax:x.set_xscale('log');x.set_yscale('log');x.grid(which='both',alpha=.3);x.axvline(t[-1],color='k',lw=.8);x.legend(fontsize=8)
fig.tight_layout();fig.savefig('f02_volume_time.png',dpi=130,facecolor='white');print('t0 Gyr',t[-1])
