# rev7 regeneration of rev6 Fig 23 (parameter degeneracies). rev6 used the R=4 Gpc stage-A fit (deg.py); this recomputes the
# same two slices on the rev7 p8 field fit at D=34.96 Gpc, R=R_ls=0.18957 Gpc: (left) window centre chi_w vs amplitude lgA,
# (right) slope n vs lgA, all other 13 parameters held at best fit. Statistic = the p8 fit measure (A.chi2_fit, l=2-2500).
import os,sys,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from multiprocessing import Pool
H=os.path.dirname(os.path.abspath(__file__));sys.path.insert(0,H);from _p8model import P8,XB,J
R=J['R_Gpc'];N=int(os.environ.get('NG','21'))
cw=np.unique(np.r_[np.linspace(0.12,R,N),XB[6]]);lA=XB[0]+np.linspace(-1.2,1.2,N);ns=XB[1]+np.linspace(-0.8,0.8,N)
def f1(t):x=XB.copy();x[6],x[0]=t;return P8.obj(x)
def f2(t):x=XB.copy();x[1],x[0]=t;return P8.obj(x)
if __name__=='__main__':
    with Pool(int(os.environ.get('NP','7'))) as p:
        S1=np.array(p.map(f1,[(c,a) for a in lA for c in cw])).reshape(N,len(cw));S2=np.array(p.map(f2,[(n,a) for a in lA for n in ns])).reshape(N,N)
    c0=P8.obj(XB);np.savez(H+'/R23_degeneracy.npz',cw=cw,lA=lA,ns=ns,S1=S1,S2=S2,c0=c0)
    fig,ax=plt.subplots(1,2,figsize=(12,4.8))
    for a,X,S,xl,xb in ((ax[0],cw,S1,'window centre $\\chi_w$ (Gpc)',XB[6]),(ax[1],ns,S2,'$P_\\delta$ slope $n$',XB[1])):
        S=np.where(S>1e14,np.nan,S);cs=a.contourf(X,lA,np.log10(S),30,cmap='viridis');plt.colorbar(cs,ax=a,label='log$_{10}$ fit $\\chi^2$')
        a.contour(X,lA,S-c0,[2.3,6.18,11.8],colors='w',linestyles=['-','--',':']);a.plot(xb,XB[0],'r*',ms=13,label=f'p8 best fit, $\\chi^2$={c0:.1f}')
        a.set_xlabel(xl);a.set_ylabel('lgA ($P_\\delta$ amplitude)');a.legend(loc='lower center',fontsize=8)
    ax[0].set_title('radial–amplitude slice (Δχ² = 2.3/6.2/11.8)',fontsize=10);ax[1].set_title('slope–amplitude slice',fontsize=10)
    fig.suptitle('Degeneracies of the rev7 field fit (D=34.96 Gpc, R=0.18957 Gpc; other 13 parameters fixed)',fontsize=10)
    fig.tight_layout();fig.savefig(H+'/R23_degeneracy.png',dpi=130);print('c0',c0,'S1min',np.nanmin(S1),'S2min',np.nanmin(S2))
