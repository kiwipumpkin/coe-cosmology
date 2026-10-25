import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
lon=np.linspace(-np.pi,np.pi,360); lat=np.linspace(-np.pi/2,np.pi/2,180)
L,B=np.meshgrid(lon,lat)
def Q(L,B,rot=0.0): return 0.5*(3*np.sin(B)**2-1)+0.6*np.cos(B)**2*np.cos(2*(L-rot))
A=1e-3
Q1=A*Q(L,B); Q2=-A*0.995*Q(L,B,rot=np.radians(0.3))
S=Q1+Q2
w=np.cos(B); rms=lambda f: np.sqrt((f**2*w).sum()/w.sum())
r1,r2,rs=rms(Q1),rms(Q2),rms(S)
fig=plt.figure(figsize=(18,10))
for i,(f,t,v) in enumerate([(Q1,f'Q₁   rms = {r1:.2e}',A),(Q2,f'Q₂ ≈ −0.995·Q₁ (rotated 0.3°)   rms = {r2:.2e}',A),(S,f'Q₁+Q₂, same scale   rms = {rs:.1e}',A),(S,f'Q₁+Q₂, zoomed to ±2×10⁻⁵',2e-5)]):
    ax=fig.add_subplot(2,2,i+1,projection='mollweide')
    im=ax.pcolormesh(L,B,f,cmap='RdBu_r',vmin=-v,vmax=v,shading='auto')
    ax.set_title(t,fontsize=16); ax.set_xticklabels([]); ax.set_yticklabels([])
    cb=fig.colorbar(im,ax=ax,orientation='horizontal',fraction=0.05,pad=0.04); cb.set_label('ΔT/T',fontsize=13); cb.ax.tick_params(labelsize=11)
fig.suptitle(r'Two nearly opposite quadrupoles: $\Delta T/T=Q_1+Q_2$, each $\sim10^{-3}$, sum $\sim10^{-5}$ (illustration, not the real sky)',fontsize=17)
plt.tight_layout(rect=[0,0,1,0.95]); plt.savefig('near_cancel_quadrupoles.png',dpi=100)
print(r1,r2,rs,np.abs(S).max())
