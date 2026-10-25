import json,os,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
os.chdir('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/constrained')
Pl=np.loadtxt('../../../reproduce_cmb/planckTT.txt');Pl=Pl[(Pl[:,0]>=2)&(Pl[:,0]<=30)]
th=np.loadtxt('/tmp/th.txt');th=th[(th[:,0]>=2)&(th[:,0]<=30)]
S=(1,.5,.25,.1,.05,.02,.01,.005,.002,.001)
rows=[('0.23',3,'red_surface'),('1',7,'R1'),('2',6,'R2'),('3',5,'R3'),('4',4,'R4'),('14.3',3,'paper_R14.3')]
rows=[r for r in rows if os.path.exists(f'red_web{r[1]}_results.json')]
fig,ax=plt.subplots(len(rows),1,figsize=(10,3.6*len(rows)),sharex=True);l=np.arange(31);cm=plt.cm.viridis
md=["| distance (Gpc) | best strength | CMB χ² | D2 | D3 | D4 | D5 | D6 | max D_ℓ ℓ7–30 |","|---|---|---|---|---|---|---|---|---|"];grid=[]
for a,(name,n,tag) in zip(ax,rows):
  d=json.load(open(f'red_web{n}_results.json'))
  a.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='ko',ms=3,label='Planck')
  a.plot(th[:,0],th[:,1],'k-',lw=1.5,label='ΛCDM χ²=20.1')
  sm=d['smooth'];a.plot(l[2:],np.array(sm['Dl'])[2:],'r--',lw=1.5,label=f"smooth red χ²={sm['chi2']:.1f}")
  cs=[]
  for i,s in enumerate(S):
    r=d[f'{tag}_{s}'];cs.append(r['chi2']);a.plot(l[2:],np.maximum(np.array(r['Dl'])[2:],1e-3),color=cm(i/9),lw=1,label=f'{s*100:g}% χ²={r["chi2"]:.4g}')
  grid.append((name,cs));b=int(np.argmin(cs));D=d[f'{tag}_{S[b]}']['Dl']
  md.append(f"| {name} | {S[b]*100:g}% | {cs[b]:.4g} | "+" | ".join(f"{D[k]:.4g}" for k in range(2,7))+f" | {max(D[7:31]):.4g} |")
  a.set_yscale('log');a.set_ylim(1e-3,1e7);a.set_title(f'red + web, kernel at {name} Gpc');a.set_ylabel('D_ℓ (μK²)');a.legend(fontsize=6.5,ncol=3,loc='upper right')
ax[-1].set_xlim(1.5,30.5);ax[-1].set_xlabel('ℓ');fig.tight_layout();fig.savefig('red_web_stack.png',dpi=110)
md+=["","Full CMB χ² (ℓ 2–30) grid","","| distance (Gpc) | "+" | ".join(f'{s*100:g}%' for s in S)+" |","|---|"+"---|"*len(S)]
md+=[f"| {n} | "+" | ".join(f'{c:.4g}' for c in cs)+" |" for n,cs in grid]
open('red_web_stack_table.md','w').write("\n".join(md)+"\n");print("\n".join(md))
