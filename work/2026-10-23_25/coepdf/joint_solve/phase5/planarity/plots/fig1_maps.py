import sys;sys.path.insert(0,'.');from common import *
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
def P(l,b,**k):hp.projplot(l,b,lonlat=True,**k)
fig=plt.figure(figsize=(16,13))
for i,k in enumerate(ORDER):
    a=ALM[k]
    for j,(ls,t) in enumerate([([2],'l=2'),([3],'l=3'),([2,3],'l=2+3')]):
        f=np.zeros(hp.Alm.getsize(LM),complex)
        for l in ls:
            for m in range(l+1):f[hp.Alm.getidx(LM,l,m)]=a[hp.Alm.getidx(LM,l,m)]
        mp=hp.alm2map(f,64,lmax=LM)*T0;v=np.abs(mp).max()
        hp.mollview(mp,fig=fig.number,sub=(4,3,3*i+j+1),cmap='RdBu_r',min=-v,max=v,unit='muK',title=f'{LAB[k]}  {t}',notext=True)
        hp.graticule(dpar=30,dmer=60,alpha=.3)
        P(*COE,marker='*',ms=13,color='k',ls='');P(0,15,marker='*',ms=13,mfc='none',mec='k',ls='')
        P(260,60,marker='X',ms=11,color='lime',mec='k',ls='');P(80,-60,marker='X',ms=11,mfc='none',mec='g',ls='')
        # AoE-normal great circle
        n=np.array(hp.ang2vec(260,60,lonlat=True));u=np.cross(n,[0,0,1.]);u/=np.linalg.norm(u);w=np.cross(n,u);ph=np.linspace(0,2*np.pi,400)
        c=np.outer(np.cos(ph),u)+np.outer(np.sin(ph),w);lo,la=hp.vec2ang(c,lonlat=True);P(lo,la,color='g',ls='',marker='.',ms=1.5)
        for l,mk,cl in [(2,'o','m'),(3,'s','orange')]:
            ev,e=eig(a,l)
            if k=='CoE_bg_only':
                x=lb(e[:,0]);P(*x,marker='D',ms=8,mfc='none',mec=cl,ls='',mew=2)
            else:
                x=lb(e[:,-1]);P(*x,marker=mk,ms=8,mfc='none',mec=cl,mew=2,ls='');P(x[0]+180,-x[1],marker=mk,ms=8,mfc='none',mec=cl,mew=1,ls='')
        if k=='CoE_bg_only':plt.text(0,-1.25,'axisymmetric: diamonds = min-<Lz^2> (symmetry) axis',ha='center',fontsize=8)
from matplotlib.lines import Line2D as L2
h=[L2([],[],marker='*',color='k',ls='',ms=12,label='CoE (180,-15)'),L2([],[],marker='*',mfc='none',mec='k',ls='',ms=12,label='anti-CoE (0,+15)'),
 L2([],[],marker='X',color='lime',mec='k',ls='',ms=10,label='AoE (260,+60)'),L2([],[],marker='X',mfc='none',mec='g',ls='',ms=10,label='AoE antipode'),
 L2([],[],color='g',ls=':',label='AoE-normal great circle (preferred plane)'),L2([],[],marker='o',mfc='none',mec='m',ls='',mew=2,label='l=2 max-<Lz^2> axis (thick) & antipode'),
 L2([],[],marker='s',mfc='none',mec='orange',ls='',mew=2,label='l=3 max-<Lz^2> axis & antipode'),L2([],[],marker='D',mfc='none',mec='m',ls='',mew=2,label='CoE bg: l=2/3 min-axis (m/orange)')]
fig.legend(handles=h,loc='lower center',ncol=4,fontsize=9);fig.suptitle('Galactic Mollweide, l increases to the left; healpy default')
plt.savefig('fig1_lowl_maps.png',dpi=110,bbox_inches='tight')
