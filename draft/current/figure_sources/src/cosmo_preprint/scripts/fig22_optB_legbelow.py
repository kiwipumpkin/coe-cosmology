import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
BLUE,ORANGE,GREEN='#1f6fc4','#e07b10','#2a9a2a'
def curves(u):
    """u in [0,1] = schematic redshift coordinate. C = G*(1-f), residual = G*f = |G-C| exactly."""
    s=lambda c,w:1/(1+np.exp(-(u-c)/w))
    n=lambda c,w:(s(c,w)-1/(1+np.exp(c/w)))/(1/(1+np.exp(-(1-c)/w))-1/(1+np.exp(c/w)))
    G=10**(-0.15-1.85*n(0.45,0.15))          # geometry: ~0.7 -> 1e-2
    f=(1-1e-3)*10**(-3*n(0.55,0.12))          # residual fraction: ~1 -> 1e-3
    C=G*(1-f); R=np.abs(G-C)
    assert np.all(C<G) and np.all(R<=np.maximum(G,C)+1e-15)
    return G,C,R
LG,LC,LR='Geometry $G(z)$','Correction $C(z)$\n(SR + mass)','Residual\n$|G-C|$'
def plot3(ax,x,G,C,R):
    ax.plot(x,G,'-',c=BLUE,lw=2,label=LG); ax.plot(x,C,'--',c=ORANGE,lw=2,label=LC)
    ax.plot(x,R,'-.',c=GREEN,lw=2,label=LR)
def floor(ax,label=True):
    ax.axhline(1e-5,ls=(0,(4,3)),c='0.35',lw=1.1)
    if label: ax.text(0.005,1.3e-5,'CMB floor $\\sim10^{-5}$',transform=ax.get_yaxis_transform(),fontsize=9.5,color='0.3',va='bottom')
z=np.logspace(-2,np.log10(1100),800); x=np.log10(z); X=lambda zz:np.log10(zz)
u=(x-x.min())/(x.max()-x.min()); G,C,R=curves(u)
fig,ax=plt.subplots(figsize=(8.4,2.95),dpi=220)
edges=[x.min(),X(0.3),X(3),X(30),x.max()]; cols=['#dceaf9','#ece3f6','#fde6d3','#dff0dc']
labs=[('Low z: Migkas\nanisotropy persists','#1a5aa8'),('Intermediate:\nmatter dipole','#5b3a8c'),
      ('Transition:\npartial cancel.','#a8441a'),('High z: strong\ncancel. → floor','#1f7a1f')]
for i in range(4):
    ax.axvspan(edges[i],edges[i+1],color=cols[i],alpha=.85,lw=0)
    ax.text((edges[i]+edges[i+1])/2,12,labs[i][0],ha='center',va='center',fontsize=7.8,fontweight='bold',color=labs[i][1])
plot3(ax,x,G,C,R); floor(ax,False); ax.text(X(3.3),1.3e-5,r'CMB floor $\sim10^{-5}$',fontsize=9,color='0.3',va='bottom')
for zz,t in[(0.3,'Migkas\nclusters\n$z\\lesssim0.3$'),(1,'CatWISE\nquasars\n$z\\sim1$'),(1100,'CMB\n$z\\sim1100$')]:
    xx=X(zz); ax.axvline(xx-(0.01 if zz>100 else 0),ymin=0,ymax=0.04,c='k',lw=2)
    ax.text(xx+(-0.04 if zz>100 else 0.04),(2e-4 if zz>100 else 1.8e-5),t,ha='right' if zz>100 else 'left',va='bottom',fontsize=8)
    ax.axvline(xx,c='k',lw=.6,ls=':',ymax=0.68)
ax.axvspan(x.min(),X(0.3),ymin=0,ymax=0.02,color='k',alpha=.5)
ax.set_yscale('log'); ax.set_ylim(1e-6,600); ax.set_xlim(x.min(),x.max())
ax.set_ylabel('Rel. amplitude\n(log, illustr.)',fontsize=9)
tz=[0.01,0.1,0.3,1,3,10,100,1100]; ax.set_xticks([X(t) for t in tz]); ax.set_xticklabels([str(t) for t in tz],fontsize=9.5)
ax.set_xlabel('Redshift $z$ (log axis, schematic)',fontsize=9); ax.tick_params(axis='y',labelsize=9.5)
ax.set_title('Relative amplitudes vs. redshift (schematic, not to scale)',fontsize=10,fontweight='bold')
ax.legend(loc='upper center',bbox_to_anchor=(0.5,-0.22),ncol=5,fontsize=8,frameon=False,handlelength=2.2)
fig.tight_layout(); fig.savefig('/workspace/cosmo_preprint/figs/Fig_SR_z_optB_legbelow.png',bbox_inches='tight',pad_inches=0.05)
