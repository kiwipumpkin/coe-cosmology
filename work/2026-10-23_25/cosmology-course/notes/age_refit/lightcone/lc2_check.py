import json,pickle,sys,numpy as np; from lc2 import *
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
r=json.load(open('res_exp.json')); L,tn=make2('exp',r['x'],Ri=1e-3)
ts=tstar(L,tn); o=obs(L,tn); lam=o['Hbar']*F.C/(70/F.KMSMPC); X0=np.array([0.]); RO=lambda t: L.R(t,X0)[0]
print('tn',tn,'t*',ts,'U/c',o['U'],flush=True)
fin=emitters_circle(L,tn,ts,np.arange(0,360,30.)); pickle.dump(dict(fin=fin,ts=ts,tn=tn),open('lc2_check.pkl','wb'))
fig,ax=plt.subplots(figsize=(18,13)); ax.set_aspect('equal')
tt=np.geomspace(ts,tn,200); ax.plot([RO(t)*lam for t in tt],np.zeros(200),'k-',lw=3,label='observer path O(t)')
cols=['red','orange','gold','green','blue','indigo','violet']; th=np.linspace(0,2*np.pi,300)
for k,t in enumerate(np.linspace(ts,tn,8)[1:]):
    ax.plot((RO(t)+(t-ts)*np.cos(th))*lam,(t-ts)*np.sin(th)*lam,c=cols[k],lw=2.5,label=f'C{k+1}: radius c(t−t*) about O(t), t={t/tn:.2f} tn'); ax.plot(RO(t)*lam,0,'*',c=cols[k],ms=20,mec='k')
for e in fin:
    p=e['pts']; Rr=L.R(p['t'],p['x']); X=Rr*np.cos(p['phi'])*lam; Y=Rr*np.sin(p['phi'])*lam
    ax.plot(X,Y,'-',c='gray' if e['status']=='emitter' else 'red',lw=2)
    if 'R_e' in e:
        ax.plot(X[-1],Y[-1],'o',c='C0' if e['status']=='emitter' else 'red',ms=12)
        ax.plot((e['R_o_te']+(e['t_e']-ts)*np.cos(th))*lam,(e['t_e']-ts)*np.sin(th)*lam,':',c='gray',lw=1)
    else: ax.plot(X[-1],Y[-1],'kx',ms=14,mew=3)
ax.plot(0,0,'rx',ms=25,mew=4); ax.text(0.2,-0.6,'CoE',fontsize=22,color='r')
nb=sum(e['status']=='beyond_CoE' for e in fin); ne=sum(e['status']=='emitter' for e in fin); nst=sum(e['status'].startswith('stopped') for e in fin)
ax.set_title(f'exact-GR rays (exp best fit, U={o["U"]:.2f}c): today\'s 12 directions → emitter on circle C(t_e)\nemitters {ne}, beyond CoE {nb}, integrator-stopped {nst}  (dotted = each emitter\'s own circle)',fontsize=21)
ax.set_xlabel('Gpc along CoE→us axis (areal R cos φ)',fontsize=20); ax.set_ylabel('Gpc',fontsize=20); ax.tick_params(labelsize=18); ax.legend(fontsize=12,loc='lower left'); ax.grid()
fig.tight_layout(); fig.savefig('lc2_check.png',dpi=80)
print('FIN',[(e['ell'],e['status'],round(e.get('t_e',np.nan)/tn,4),round(e.get('z',np.nan),4),round(e.get('T_obs',np.nan),3),e.get('plasma'),e.get('no_matter')) for e in fin],flush=True)
