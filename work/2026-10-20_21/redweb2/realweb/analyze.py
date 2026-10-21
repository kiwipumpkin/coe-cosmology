import json,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb';Pl=np.loadtxt(B+'/planckTT.txt');Pl=Pl[(Pl[:,0]>=2)&(Pl[:,0]<=30)]
T0=2.7255e6;ll=np.arange(31);l=np.arange(2,31)
def Dl(x):return ll*(ll+1)*hp.anafast(np.asarray(hp.remove_dipole(x)),lmax=30)*T0**2/(2*np.pi)
def chi(D):
  c=0
  for r in Pl:
    L=int(r[0]);e=r[3] if D[L]>r[1] else r[2];c+=((D[L]-r[1])/e)**2
  return float(c)
r2=lambda D:float(np.corrcoef(D[2:31],Pl[:,1])[0,1]**2)
O='/workspace/redweb2/realweb/';z=np.load(O+'groups.npz');S=z['S'].sum(1);NG=len(S);RS=[0.5,0.9,2.0,4.0,8.0,14.0]
sm=np.load('/workspace/redweb2/sm_red.npy');E=[.01,.02,.05,.1,.2,.5,1.0];rng=np.random.default_rng(42)
res=dict(n_groups=NG,smooth=dict(chi2=chi(Dl(sm)),r2=r2(Dl(sm)),Dl=Dl(sm).tolist()),R={})
def ent(w):D=Dl(sm+w);return dict(chi2=chi(D),r2=r2(D),Dl=D.tolist())
for iR,R in enumerate(RS):
  Sg=S[:,iR];tot=Sg.sum(0);d=dict(web_rms=float(tot.std()),smooth_rms=float(sm.std()),a={},b={},c={})
  for e in E:
    k=max(1,int(round(e*NG)));d['b'][str(e)]=ent(e*tot)
    d['a'][str(e)]=[dict(ent(Sg[rng.choice(NG,k,replace=False)].sum(0)*NG/NG),n_groups=k) for _ in range(5)] if e<1 else [ent(tot)]
    cl=[]
    for _ in range(5):
      keep=rng.uniform(e,1);kk=max(1,int(round(keep*NG)));sc=e/(kk/NG);cl.append(dict(ent(sc*Sg[rng.choice(NG,kk,replace=False)].sum(0)),keep=kk/NG,scale=sc))
    d['c'][str(e)]=cl
  res['R'][str(R)]=d;print(R,'rms web/smooth %.3g'%(tot.std()/sm.std()),{m:{e:round(min(x['chi2'] for x in (d[m][e] if isinstance(d[m][e],list) else [d[m][e]])),1) for e in d[m]} for m in 'abc'},flush=True)
json.dump(res,open(O+'realweb_results.json','w'))
# plots
cols=plt.cm.tab10.colors
for R in map(str,RS):
  d=res['R'][R];cur=[('smooth red (no web)',np.array(res['smooth']['Dl']),res['smooth']['chi2'],'k','--')]
  for i,e in enumerate(E):
    cur.append((f'(b) scale {e*100:g}%',np.array(d['b'][str(e)]['Dl']),d['b'][str(e)]['chi2'],cols[i],'-'))
    for m,ls in (('a',':'),('c','-.')):
      L=d[m][str(e)];bb=min(L,key=lambda x:x['chi2']);cur.append((f'({m}) {e*100:g}% best of {len(L)}',np.array(bb['Dl']),bb['chi2'],cols[i],ls))
  for kind in ('log','lin'):
    plt.figure(figsize=(13,8))
    for lab,D,c2,c,ls in cur:
      vis=kind=='log' or np.any((D[2:31]>0)&(D[2:31]<3000))
      plt.plot(l,D[2:31],color=c,ls=ls,lw=1.6 if ls=='-' else 1.2,marker='.' ,label=(f'{lab}: χ²={c2:.4g}, r²={r2(D):.2f}' if vis else None))
    plt.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='r',ms=4,capsize=2,label='Planck TT')
    if kind=='log':plt.yscale('log');plt.ylim(bottom=1e-3)
    else:plt.ylim(0,3000)
    plt.title(f'Realistic web, sampling R={R} Gpc: (a) random removal, (b) force scale, (c) random mix; ΛCDM χ²=20.1');plt.xlabel('ℓ');plt.ylabel('D_ℓ (μK²)')
    plt.legend(fontsize=6,ncol=3,loc='lower center',framealpha=0.85);plt.tight_layout();plt.savefig(O+f'realweb_R{R}_{kind}.png',dpi=100);plt.close()
fig,ax=plt.subplots(1,3,figsize=(18,5))
for a,m in zip(ax,'abc'):
  H=np.array([[min(x['chi2'] for x in (res['R'][str(R)][m][str(e)] if isinstance(res['R'][str(R)][m][str(e)],list) else [res['R'][str(R)][m][str(e)]])) for e in E] for R in RS])
  im=a.imshow(np.log10(H),cmap='viridis_r',aspect='auto');a.set_xticks(range(len(E)),[f'{e*100:g}%' for e in E]);a.set_yticks(range(len(RS)),[f'{R}' for R in RS])
  for i in range(len(RS)):
    for j in range(len(E)):a.text(j,i,f'{H[i,j]:.3g}',ha='center',va='center',fontsize=7,color='w' if np.log10(H[i,j])>np.log10(H).mean() else 'k')
  a.set_xlabel('effect size e');a.set_ylabel('sampling R (Gpc)');a.set_title({'a':'(a) random removal, best of 5','b':'(b) force scale','c':'(c) random mix, best of 5'}[m]+' — χ²');plt.colorbar(im,ax=a,label='log10 χ²')
plt.tight_layout();plt.savefig(O+'realweb_chi2_heatmap.png',dpi=110)
