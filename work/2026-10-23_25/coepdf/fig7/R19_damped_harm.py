# rev7 regeneration of rev6 Fig 19 (damped-harmonic fits to the TT residual). rev6 fitted the residual of the R=4 Gpc acoustic fit.
# Here eps_l=(D^Planck_l - D^CoE_l)/sigma_l uses the rev7 ensemble-mean spectrum (p6_p8: D=34.96 Gpc LTB background, p8 field,
# p5_c_p8 bodies). Models (rev6 Eqs. 22-23): undamped sum_j a_j sin(2 pi j l/DL+phi_j); damped exp(-l/L) x same. Linear least
# squares in (a cos phi, a sin phi) at fixed (DL,L); DL and L scanned on a grid.
import os,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
H=os.path.dirname(os.path.abspath(__file__));Z=np.load('/workspace/coepdf/joint_solve/phase2/p6_p8.npz')
Pf=Z['Pf'];l=Pf[:,0];sig=0.5*(Pf[:,2]+Pf[:,3]);D=Z['D'][l.astype(int)];eps=(Pf[:,1]-D)/sig
DLs=np.arange(150,1601,10.);Ls=np.r_[np.arange(100,3001,50.),1e9]
def fit(n,DL,L):
    X=np.column_stack([f(2*np.pi*j*l/DL) for j in range(1,n+1) for f in (np.sin,np.cos)])*np.exp(-l/L)[:,None]
    b=np.linalg.lstsq(X,eps,rcond=None)[0];m=X@b;return 1-np.sum((eps-m)**2)/np.sum((eps-eps.mean())**2),m,np.hypot(b[0::2],b[1::2])
res={}
for n in range(1,10):
    for tag,LL in (('damped',Ls[:-1]),('undamped',[1e9])):
        best=max(((fit(n,a,L)[0],a,L) for a in DLs for L in LL));res[f'{tag}{n}']=dict(r2=best[0],DL=best[1],L=best[2] if best[2]<1e8 else None)
json.dump(res,open(H+'/R19_damped_harm.json','w'),indent=1)
fig=plt.figure(figsize=(13,12));g=fig.add_gridspec(4,1,height_ratios=[1,1,1,0.9])
for i,n in enumerate((7,8,9)):
    r=res[f'damped{n}'];r2,m,a=fit(n,r['DL'],r['L'] or 1e9);ax=fig.add_subplot(g[i])
    ax.plot(l,eps,'.',color='grey',ms=2,label='ε_ℓ = (Planck − CoE)/σ (rev7 mean)');ax.plot(l,m,'r',lw=1.2,label=f'damped, {n} harmonics: Δℓ={r["DL"]:.0f}, L={r["L"]:.0f}, r²={r2:.3f}')
    ax.axhline(0,color='k',lw=.5);ax.set_ylabel('ε_ℓ (σ)');ax.legend(fontsize=8,loc='upper right');ax.set_xlim(2,2500)
ax.set_xlabel('ℓ');ax=fig.add_subplot(g[3]);nn=np.arange(1,10)
ax.plot(nn,[res[f'damped{n}']['r2'] for n in nn],'ro-',label='damped');ax.plot(nn,[res[f'undamped{n}']['r2'] for n in nn],'bs--',label='undamped')
ax.set_xlabel('number of harmonics');ax.set_ylabel('Pearson r²');ax.legend();ax.grid(alpha=.3)
fig.suptitle('Damped-harmonic decomposition of the rev7 TT residual (D=34.96 Gpc, R=0.18957 Gpc)',fontsize=11)
fig.tight_layout();fig.savefig(H+'/R19_damped_harm.png',dpi=110);print(json.dumps(res))
