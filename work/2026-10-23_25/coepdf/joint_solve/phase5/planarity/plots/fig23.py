import sys,json;sys.path.insert(0,'.');from common import *
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
rng=np.random.default_rng(0);NS=10000;st={l:[] for l in range(2,LM+1)};sz={l:[] for l in range(2,LM+1)}
for i in range(NS):
    for l in range(2,LM+1):
        p=rng.normal(size=l+1)+1j*rng.normal(size=l+1);p[0]=p[0].real*np.sqrt(2)
        v=np.r_[[(-1)**m*np.conj(p[m]) for m in range(l,0,-1)],p];q=np.abs(v)**2;st[l].append((q[0]+q[-1])/q.sum());sz[l].append(q[l]/q.sum())
ls=range(2,LM+1);out={}
for k in ORDER:out[k]={l:dict(zip(['t_AoE','zon_AoE'],tl_about(ALM[k],AOE,l)))|dict(zip(['t_CoE','zon_CoE'],tl_about(ALM[k],COE,l))) for l in ls}
fig,ax=plt.subplots(1,2,figsize=(13,5))
for a,(key,sim,yl) in zip(ax,[('t_AoE',st,'planarity t_l about AoE'),('zon_AoE',sz,'zonal fraction |a_l0|^2/sum about AoE')]):
    q=np.array([np.percentile(sim[l],[2.5,16,84,97.5]) for l in ls]);x=np.array(list(ls))
    a.fill_between(x,q[:,0],q[:,3],color='0.85',label='isotropic sims 95%');a.fill_between(x,q[:,1],q[:,2],color='0.65',label='68%')
    for i,k in enumerate(ORDER):a.plot(x+(i-1.5)*.06,[out[k][l][key] for l in ls],'o-',label=LAB[k])
    a.set_xticks(x);a.set_xlabel('l');a.set_ylabel(yl);a.legend(fontsize=8)
fig.suptitle(f'{NS} isotropic Gaussian sims (frame-independent)');plt.savefig('fig2_planarity_zonal.png',dpi=120,bbox_inches='tight')
fig,ax=plt.subplots(1,2,figsize=(12,5.5))
for a,key,nm in [(ax[0],'t','planarity t_l'),(ax[1],'zon','zonal fraction')]:
    for i,k in enumerate(ORDER):
        x=[out[k][l][key+'_CoE'] for l in ls];y=[out[k][l][key+'_AoE'] for l in ls];a.plot(x,y,'o-',label=LAB[k],alpha=.8)
        for l,xx,yy in zip(ls,x,y):a.annotate(str(l),(xx,yy),fontsize=8)
    a.plot([0,1],[0,1],'k:',lw=.8);a.set_xlabel(nm+' about CoE axis (180,-15)');a.set_ylabel(nm+' about AoE (260,60)');a.legend(fontsize=8);a.set_xlim(0,1);a.set_ylim(0,1)
fig.suptitle('labels = l');plt.savefig('fig3_AoE_vs_CoE.png',dpi=120,bbox_inches='tight')
json.dump({k:{str(l):v for l,v in d.items()} for k,d in out.items()},open('fig23_values.json','w'),indent=1)
