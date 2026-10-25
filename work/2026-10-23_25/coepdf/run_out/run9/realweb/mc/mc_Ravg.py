import numpy as np,json
from multiprocessing import Pool
exec(open('/workspace/redweb2/realweb/analyze.py').read().split('O=')[0])
S=np.load('mc_sky.npy',mmap_mode='r');sm=np.load('/workspace/redweb2/sm_red.npy');RS=[0.5,0.9,2.0,4.0,8.0,14.0]
E=[1,.5,.2,.1,.05,.02,.01,.005,.002,.0005];pm=np.median(Pl[2:,1]);G=json.load(open('samples/bracket_guesses.json'))
def f(a):
  k,iR=a;sk=np.array(S[k][iR]);st=[];Ds=[]
  for e in E:
    D=Dl(sm+e*sk);y=D[4:31];Ds.append(D);st.append((np.median(y),abs(np.polyfit(l[2:],np.log(y),1)[0])<0.05))
  i=[i for i in range(9) if st[i][1] and st[i+1][1] and st[i][0]>=pm>=st[i+1][0]][0]
  return (Ds[i]+Ds[i+1])/2
if __name__=='__main__':
  z=np.load('mc_mean_passing.npz');RR=[4.0,8.0,14.0]
  ks={R:[k for k,o in enumerate(G) if o[str(R)]] for R in RR}
  jobs=[(k,RS.index(R)) for R in RR for k in ks[R]]
  with Pool() as p: br=p.map(f,jobs,chunksize=8)
  A={};B={};n=0
  for R in RR:
    for j,k in enumerate(ks[R]):
      A.setdefault(k,[]).append(br[n]);B.setdefault(k,[]).append(z[f'D{R}'][j]);n+=1
  U=sorted(A);a=np.array([np.mean(A[k],0) for k in U]);b=np.array([np.mean(B[k],0) for k in U])
  np.savez('mc_Ravg.npz',U=U,a=a,b=b,nR=[len(A[k]) for k in U])
  t=np.loadtxt('/workspace/joint/lcdm.txt');Lt=np.interp(l,t[:,0],t[:,1])
  fig,ax=plt.subplots(2,1,figsize=(11,11),sharex=True)
  for X,D,tt in zip(ax,[a,b],['(a) mean of bracketing ladder pair, averaged over passing R','(b) interpolated E* spectrum, averaged over passing R']):
    m=D.mean(0);lo,hi=D.min(0),D.max(0);p16,p84=np.percentile(D,[16,84],0);cl=[L for L in range(2,31) if hi[L]>2500]
    print(tt,f'N={len(D)} chi2={chi(m):.1f} r2={r2(m):.3f}');print('  mean',np.round(m[2:6]),'p16',np.round(p16[2:6]),'p84',np.round(p84[2:6]),'min',np.round(lo[2:6]),'max',np.round(hi[2:6]),'clip',cl)
    X.fill_between(l,lo[2:31],hi[2:31],color='lightblue',alpha=.6,label=f'min–max of {len(D)} universes')
    X.fill_between(l,p16[2:31],p84[2:31],color='tab:blue',alpha=.35,label='16–84 percentile')
    X.plot(l,m[2:31],color='tab:blue',lw=2,label=f'mean: χ²={chi(m):.1f}, r²={r2(m):.2f}')
    X.plot(l,Lt,'k-',lw=1.6,label='ΛCDM: χ²=20.1')
    X.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='k',mfc='w',ms=5,capsize=2,label='Planck TT')
    X.set_ylim(0,2500);X.set_ylabel('D_ℓ [μK²]');X.set_title(tt+(f'  (max clipped at ℓ={",".join(map(str,cl))})' if cl else ''),fontsize=10)
    X.legend(loc='lower center',ncol=3,fontsize=9)
  ax[-1].set_xlabel('ℓ');plt.tight_layout();plt.savefig('mc_Ravg_lin.png',dpi=120)
  nR=np.array([len(A[k]) for k in U]);print('nR counts',{c:int((nR==c).sum()) for c in (1,2,3)})
