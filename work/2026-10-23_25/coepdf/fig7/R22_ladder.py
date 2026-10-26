# rev7 regeneration of rev6 Fig 22 (same source parameters at R=4,5,6 Gpc). At the rev7 background R is not free: R=R_ls=eta_rec*D
# =0.18957 Gpc. This recomputes the field spectrum with the p8 best-fit source parameters at R = 0.8, 1.0, 1.2 x R_ls (window kept
# at the same depth R-chi_w below the shell), showing how the peak positions move with shell radius. Field term + Dfix, CORR applied.
import os,sys,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
H=os.path.dirname(os.path.abspath(__file__));sys.path.insert(0,H);from _p8model import P8,XB,J,M,A
R0=J['R_Gpc'];s=R0/4;p0=P8.P(XB);dep=R0-p0['chiw'];ll=np.arange(M.LMAX+1)
def pk(D):
    from scipy.signal import find_peaks
    i,_=find_peaks(D[100:1800],prominence=0.03*D[100:1800].max());return [int(100+j) for j in i][:5]
out={};fig,axs=plt.subplots(1,3,figsize=(18,5.2),sharey=True)
for ax,f in zip(axs,(0.8,1.0,1.2)):
    R=f*R0;M.EX=M.AFwd(R=R,lmax=M.LX,kmax=200/s,kmin=1e-2/s,nk=500,dch=0.004*s);p=dict(p0,chiw=R-dep)
    D=M.Dfield(p)*P8.CORR+P8.Dfix;c=A.chi2_fit(D);out[f'{f:g}']=dict(R=R,chi2_fit=float(c),peaks=pk(D))
    ax.plot(ll[2:],A.Lt[2:],'k',lw=1.5,label='ΛCDM');ax.plot(A.Pf[:,0],A.Pf[:,1],'.',color='grey',ms=1.5,label='Planck 2018')
    ax.plot(ll[2:],D[2:],'r',lw=1.2,label=f'CoE field, R={R:.4f} Gpc')
    ax.set_title(f'R = {f:g}·R_ls = {R:.4f} Gpc;  fit χ²={c:.0f}\npeaks {out[f"{f:g}"]["peaks"]}',fontsize=9);ax.set_xlabel('ℓ');ax.set_xlim(2,2500);ax.legend(fontsize=7)
axs[0].set_ylabel('D_ℓ (μK²)');fig.suptitle('Same p8 source parameters at three shell radii (rev7 background D=34.96 Gpc; R_ls=0.18957 Gpc is the physical value)',fontsize=10)
fig.tight_layout();fig.savefig(H+'/R22_ladder.png',dpi=120);json.dump(out,open(H+'/R22_ladder.json','w'),indent=1);print(out)
