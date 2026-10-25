import numpy as np,json
d=json.load(open('cmbann.json')); a=d['annulus']; b,S,dop,gm,te=np.load('cmbann_keep.npy')
print('share no-gamma: stretch',np.mean(np.log(S))/np.mean(np.log(S*dop)),'motion',np.mean(np.log(dop))/np.mean(np.log(S*dop)))
Dn=(a['now'][0]+a['now'][1])/2; print('centre offset Gpc',Dn-26.8378); assert abs(Dn-26.8378)<0.01
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
fig,ax=plt.subplots(1,3,figsize=(36,12)); ph=np.linspace(0,2*np.pi,800)
for A,(ri,ro),Dp,title,scale in ((ax[0],a['at_trec'],a['D_trec'],f"At recombination, t = {d['trec_Gyr']*1e3:.2f} Myr",1e3),(ax[1],a['now'],Dn,'Now, t = 2.665 Gyr',1)):
    u='Mpc' if scale==1e3 else 'Gpc'; ri*=scale; ro*=scale; Dp*=scale
    A.fill(np.r_[ro*np.cos(ph),ri*np.cos(ph[::-1])],np.r_[ro*np.sin(ph),ri*np.sin(ph[::-1])],color='orange',alpha=.35,label='light-filled annulus (layers q_in–q_out)')
    A.plot(0,0,'*',color='r',ms=30,label='CoE + mass'); A.plot(Dp,0,'o',color='k',ms=16,label=f'us, D = {Dp:.2f} {u}')
    A.set_aspect('equal'); L=ro*1.15; A.set_xlim(-L,L); A.set_ylim(-L,L); A.tick_params(labelsize=20)
    A.set_xlabel(u,fontsize=24); A.set_title(title+f"\nedges {ri:.2f}–{ro:.2f} {u}",fontsize=26); A.legend(fontsize=18,loc='upper left')
ax[1].text(-0.95*a['now'][1]*1.1,-0.95*a['now'][1]*1.05,f"V = 4π∫r²dr\nV_now/V_emit = {a['V_ratio']:.3g}\nannulus stretch = {a['stretch']:.0f}\nR_rec = {d['R_rec_Mpc']:.1f} Mpc",fontsize=22,bbox=dict(fc='w'))
B=ax[2]; B.plot(b,3000/S,lw=4,label=f"stretch only: mean {np.mean(3000/S):.3f} K")
B.plot(b,3000/(S*dop),lw=4,label=f"stretch × Doppler: mean {np.mean(3000/(S*dop)):.3f} K")
B.axhline(2.725,color='k',ls=':',lw=2); B.set_xlim(0,360); B.tick_params(labelsize=20)
B.set_xlabel('arrival direction β (deg, 0 = away from CoE)',fontsize=24); B.set_ylabel('T_obs (K)',fontsize=24)
B.set_title('T(direction); γ undefined (v_rel ≥ c), gravity undefined',fontsize=26); B.legend(fontsize=20)
fig.tight_layout(rect=(0,0,1,0.93)); fig.savefig('cmb_annulus.png',dpi=60)
