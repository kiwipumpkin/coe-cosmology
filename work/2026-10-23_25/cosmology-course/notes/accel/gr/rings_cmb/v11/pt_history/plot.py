import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=json.load(open('pt_history.json'))['rays'];o=np.load('../clusters_voids/quad_balance/mass_family_scan/runs/soft_2_30.npz')['o']
f,ax=plt.subplots(1,3,figsize=(17,5))
c={'0':'C3','90':'C2','180':'C0'}
for k,r in d.items():
    w=np.array(r['w']);E=-np.array(r['pt']);t=np.array(r['t'])
    ax[0].semilogy(w,E,c[k],label=f"psi={k}°: 1+z={E[-1]:.2f} at w={w[-1]:.3f}");ax[0].plot(w[-1],E[-1],'o',c=c[k])
    ax[1].loglog(t,E,c[k],label=f"psi={k}°")
ax[0].axvline(1,ls=':',c='k');ax[0].text(1.0,1.5,' observer (w=1)')
ax[0].set_xlabel('CoE-distance label w  (CoE side ←  → far side)');ax[0].set_ylabel('photon energy E=-pt (observer=1)')
ax[0].set_title('Energy along each ray vs CoE distance');ax[0].legend(fontsize=8)
ax[1].invert_xaxis();ax[1].set_xlabel('cosmic time t (code units, back toward emission →)');ax[1].set_ylabel('E=-pt');ax[1].set_title('Energy vs time: all three nearly identical');ax[1].legend()
ax[2].plot(o[:,0],o[:,1],'k.-');ax[2].set_xlabel('psi (deg) of photon momentum; 0 = moving toward CoE');ax[2].set_ylabel('1+z at emission');ax[2].set_title('Net 1+z vs psi (soft p=2, frac=30)')
plt.tight_layout();plt.savefig('pt_history.png',dpi=130)
