import json,numpy as np;exec(open('jump_an.py').read().split("# per-body")[0].replace("for r in rows:print","for r in []:print"))
nv=json.load(open('thin4_results.json'))['n_voids'];full=base+W.sum(0);c1=chi(l_(full));Df=l_(full)
s=(u>=0.0170)&(u<0.0171);i0=base+cs[int((u<0.0170).sum())-1]
out=[]
for g_ in np.unique(gid[s]):
  k=np.where(gid==g_)[0];t='void (-D,+W pair)' if g_<nv else 'cluster'
  Da=l_(i0+W[k].sum(0));Dr=l_(full-W[k].sum(0));Di=l_(i0)
  out.append(dict(gid=int(g_),type=t,M_parts=M[k].tolist(),M_net=float(M[k].sum()),pos=X[k[0]].tolist(),r=float(np.linalg.norm(X[k[0]])),dsphere=float(np.linalg.norm(X[k[0]])-4),
    chi_before=chi(Di),chi_add=chi(Da),chi_full_without=chi(Dr),D30_add=float(Da[30]-Di[30]),D25_add=float(Da[25]-Di[25]),D30_without=float(Dr[30]),Dl_unit_alone=l_(W[k].sum(0)).tolist()))
for o in sorted(out,key=lambda o:-o['D30_add']):print(o['type'],'gid',o['gid'],'Mnet %.2e'%o['M_net'],np.round(o['pos'],3),'r=%.4f d4=%+.4f'%(o['r'],o['dsphere']),'chi %.1f->%.1f'%(o['chi_before'],o['chi_add']),'dD30 %.0f dD25 %.0f'%(o['D30_add'],o['D25_add']),'1.75%% without: chi %.1f D30 %.0f'%(o['chi_full_without'],o['D30_without']))
json.dump(out,open('thin4_jump_step.json','w'))
