# Mollweide of l=2+3 (Galactic) for Planck target vs v3 seed-1 model sky; marks AoE (260,+60) and best <Lz^2> axes for l=2,3 (from planarity.json).
import numpy as np,healpy as hp,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
V='/workspace/coepdf/joint_solve/phase5/web/inverse/v3/';T0=2.7255e6;d=np.load(V+'inv3_s1.npz');R=json.load(open('planarity.json'))
f=np.zeros(6);f[2:4]=1;fig=plt.figure(figsize=(14,5))
for i,(k,t) in enumerate([('Planck_target','Planck SMICA (l=2+3)'),('v3_s1_model','v3 seed-1 model sky (l=2+3)')]):
    m=hp.alm2map(hp.almxfl(hp.map2alm(d['Tobs' if i==0 else 'Resid'],lmax=5,iter=3),f),64,lmax=5)*T0
    hp.mollview(m,fig=fig.number,sub=(1,2,i+1),title=f"{t}\nt2={R[k]['2']['t_l']:.3f}, t3={R[k]['3']['t_l']:.3f}, angle(l2,l3)={R[k]['ang_l2_l3']:.1f} deg",unit='muK',cmap='RdBu_r',min=-60,max=60)
    hp.graticule(30,alpha=.3)
    for lb,mk,c,lab in [((260,60),'*','k','AoE'),(R[k]['2']['best_axis_lb'],'o','g','l=2 axis'),(R[k]['3']['best_axis_lb'],'s','m','l=3 axis')]:
        for s in (1,-1):hp.projplot((lb[0]+(0 if s>0 else 180))%360,s*lb[1],mk,color=c,ms=10,lonlat=True,mfc='none' if s<0 else c,label=lab if s>0 else None)
    plt.legend(loc='lower right',fontsize=8)
plt.savefig('planarity_mollweide.png',dpi=120,bbox_inches='tight')
