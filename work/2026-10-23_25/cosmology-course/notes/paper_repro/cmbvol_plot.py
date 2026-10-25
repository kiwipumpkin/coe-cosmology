import numpy as np,json
src=open('steps.py').read(); exec(src[:src.index('import matplotlib')])
b,S,dop,gam,grav,qe,te=np.load('cmbvol_keep.npy'); dq=2e-3
def lay(q,t): r=r_of(q,t); return r,(r_of(q+dq,t)-r_of(q-dq,t))
i=np.argmin(abs(b-0)); j=np.argmin(abs(b-180))
info=[]
for k in (i,j):
    re,dre=lay(qe[k],te[k]); rn,drn=lay(qe[k],1.0); info.append((b[k],re*DG,dre*DG,rn*DG,drn*DG,(rn**2*drn)/(re**2*dre),S[k]))
T1=3000/(S*dop); T2=3000/S; T3=3000/(S*dop*gam)
def st(T):
    A=np.c_[np.ones_like(b),np.cos(np.radians(b)),np.sin(np.radians(b))]; c,*_=np.linalg.lstsq(A,T,rcond=None); r=T-A@c
    return T.mean(),np.hypot(c[1],c[2]),r.std(),r.std()/T.mean()/1e-5
print('S*dop*gam',st(T3)); print('layers',info)
lz=np.log(S*dop*gam); print('share stretch',np.mean(np.log(S))/lz.mean(),'motion',np.mean(np.log(dop*gam))/lz.mean())
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
fig,(A,B)=plt.subplots(1,2,figsize=(30,12)); ph=np.linspace(0,2*np.pi,600)
for (bb,re,dre,rn,drn,ratio,s),c in zip(info,('tab:red','tab:blue')):
    for r,ls in ((re,'-'),(rn,'-')): A.fill_between(np.cos(ph)*0,0,0)
cols=['tab:red','tab:blue']
bb,re,dre,rn,drn,ratio,s=info[0]
A.plot(re*np.cos(ph),re*np.sin(ph),color='purple',lw=4,label=f'layer at emission: r = {re:.2f}, dr = {dre*1e3:.2f} Mpc')
A.plot(rn*np.cos(ph),rn*np.sin(ph),color='tab:red',lw=4,label=f'same layer now: r = {rn:.2f}, dr = {drn*1e3:.1f} Mpc')
A.plot(0,0,'*',color='k',ms=30); A.plot(r_of(1.0,1.0)*DG,0,'o',color='orange',ms=20,mec='k',label='us now')
A.text(-26,-24,f"V = 4π r² dr (3D spherical shell)\nV_now/V_emit = {ratio:.0f}\nstretch = (ratio)^(1/3) = {s:.2f}\nneeded for 2.725 K: 1101",fontsize=24,bbox=dict(fc='w',ec='k'))
A.set_aspect('equal'); A.set_xlim(-30,30); A.set_ylim(-30,30); A.legend(fontsize=20,loc='upper left')
A.set_xlabel('Gpc',fontsize=24); A.set_ylabel('Gpc',fontsize=24); A.tick_params(labelsize=20)
A.set_title('Emitter layer (soft-core fit, ray β=0°), true scale radii',fontsize=26)
B.plot(b,T2,lw=4,label=f'stretch only: mean {T2.mean():.0f} K')
B.plot(b,T1,lw=4,label=f'stretch × Doppler: mean {T1.mean():.0f} K')
B.plot(b,T3,lw=4,label=f'stretch × Doppler × γ: mean {T3.mean():.0f} K')
B.axhline(2.725,color='k',ls=':'); B.set_xlabel('arrival direction β (deg, 0 = away from CoE)',fontsize=24); B.set_ylabel('T_obs (K)',fontsize=24)
B.tick_params(labelsize=20); B.legend(fontsize=20); B.set_xlim(0,360)
B.set_title('T(direction); gravity term undefined (Φ/c² ≈ −2.8×10⁴)',fontsize=26)
fig.tight_layout(); fig.savefig('cmb_volume_stretch.png',dpi=65)
