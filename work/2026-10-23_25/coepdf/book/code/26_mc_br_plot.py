import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
exec(open('/workspace/redweb2/realweb/analyze.py').read().split('O=')[0])
S=np.load('mc_sky.npy',mmap_mode='r');sm=np.load('/workspace/redweb2/sm_red.npy');RS=[0.5,0.9,2.0,4.0,8.0,14.0]
t=np.loadtxt('/workspace/joint/lcdm.txt');l=np.arange(2,31);Lt=np.interp(l,t[:,0],t[:,1])
E=[1,.5,.2,.1,.05,.02,.01,.005,.002,.0005];G=np.logspace(-3.5,0.3,60)
ks=[0,int(__import__('sys').argv[1])];out={};B=json.load(open('samples/bracket_guesses.json'))
for k in ks:
  k=int(k);fig,ax=plt.subplots(6,1,figsize=(11,30));out[k]={}
  for a,iR,R in zip(ax,range(6),RS):
    sk=np.array(S[k][iR]);cs=[chi(Dl(sm+e*sk)) for e in G];eb=G[int(np.argmin(cs))];cb=min(cs)
    Db=Dl(sm+eb*sk)[2:31];y=Db[2:];slope=np.polyfit(l[2:],np.log(y),1)[0];over=Db[:3].max()/Pl[:3,1].max()
    ok=B[k][str(R)]
    out[k][R]=dict(best_e=float(eb),chi2=float(cb),slope=float(slope),lowl_over=float(over),guess='PASS' if ok else 'FAIL')
    a.plot(l,Lt,'k-',lw=2,label='ΛCDM χ²=20.1');D0=Dl(sm);a.plot(l,D0[2:31],'r--',lw=2,label=f'smooth red χ²={chi(D0):.1f}')
    for j,e in enumerate(E):
      D=Dl(sm+e*sk);a.plot(l,D[2:31],color=plt.cm.viridis(j/9),label=f'{e*100:g}% χ²={chi(D):.4g}')
    a.plot(l,Db,'m-',lw=2.5,label=f'best {eb*100:.3g}% χ²={cb:.0f}')
    a.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='k',ms=4,label='Planck')
    a.set_yscale('log');a.set_ylim(1e-3,None);a.set_ylabel('D_ℓ (μK²)')
    a.set_title(f'red + MC web #{k}, kernel at {R} Gpc — my guess: {out[k][R]["guess"]}',color='g' if ok else 'r')
    a.legend(loc='lower center',ncol=4,fontsize=7,framealpha=.85)
  ax[-1].set_xlabel('ℓ');plt.tight_layout();plt.savefig(f'samples/mc{k}_bracket.png',dpi=80);plt.close()
json.dump(out,open('samples/guesses.json','w'),indent=1)
for k in out:print(k,{R:(v['guess'],round(v['chi2']),f"{v['best_e']*100:.3g}%") for R,v in out[k].items()})
