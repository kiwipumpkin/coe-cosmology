import sys,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import fit
T0=2.7255e6;B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/'
o=np.load(B+'quad_balance/mass_family_scan/runs/soft_2_30.npz')['o'];o=o[np.argsort(o[:,0])]
# l=2 zonal value as a function of angle psi from CoE axis (as seen from us). Build from the smooth profile.
psi=np.linspace(0,180,361)
sm=np.exp(-CubicSpline(o[:,0],np.log(o[:,1]))(psi));sm=sm/sm.mean()-1
# project onto Legendre P2 to get the l=2 zonal part
mu=np.cos(np.radians(psi));P2=0.5*(3*mu**2-1)
w=np.sin(np.radians(psi))
a2=np.trapezoid(sm*P2*w,np.radians(psi))/np.trapezoid(P2*P2*w,np.radians(psi))
q=a2*P2*T0  # uK, l=2 zonal contribution vs psi
# Geometry: CoE at origin (LEFT), us at distance D to the right on +x axis.
D=58.2; Rvis=14.3
# sightline directions from us: psi measured from the CoE direction (which points -x from us)
th=np.radians(psi)
# unit direction of sightline: at psi=0 points toward CoE (-x). rotate.
dx=-np.cos(th); dy=np.sin(th)
ux,uy=D+Rvis*dx, Rvis*dy   # points on our visible bubble
fig,ax=plt.subplots(figsize=(10,8))
# color each visible-bubble point by q
import matplotlib as mpl
norm=mpl.colors.TwoSlopeNorm(vmin=q.min(),vcenter=0,vmax=q.max())
for sgn in (1,-1):
    sc=ax.scatter(ux,sgn*uy,c=q,cmap='RdBu_r',norm=norm,s=14)
cb=plt.colorbar(sc,ax=ax,label='l=2 (quadrupole) contribution  [μK]')
# CoE and us
ax.plot(0,0,'k*',ms=22);ax.annotate('CoE (origin)',(0,0),(3,3),fontsize=11,weight='bold')
ax.plot(D,0,'ko',ms=9);ax.annotate('us  (D = %.1f Gpc)'%D,(D,0),(D-6,-6),fontsize=11,weight='bold')
ax.plot([0,D],[0,0],'k--',lw=1,alpha=.6)
# visible bubble circle
ph=np.linspace(0,2*np.pi,200);ax.plot(D+Rvis*np.cos(ph),Rvis*np.sin(ph),'k:',lw=1,alpha=.5)
# arrows to coldest (toward/away CoE) and hottest (perp)
ax.annotate('COLD\n(toward CoE)',(D-Rvis,0),(D-Rvis-16,6),color='blue',fontsize=10,ha='center',arrowprops=dict(arrowstyle='->',color='blue'))
ax.annotate('COLD\n(away from CoE)',(D+Rvis,0),(D+Rvis+2,6),color='blue',fontsize=10,arrowprops=dict(arrowstyle='->',color='blue'))
ax.annotate('HOT (perp.)',(D,Rvis),(D+3,Rvis+3),color='red',fontsize=10,arrowprops=dict(arrowstyle='->',color='red'))
ax.annotate('HOT (perp.)',(D,-Rvis),(D+3,-Rvis-5),color='red',fontsize=10,arrowprops=dict(arrowstyle='->',color='red'))
ax.set_aspect('equal');ax.set_xlabel('x  [Gpc]');ax.set_ylabel('y  [Gpc]')
ax.set_title('Quadrupole on our sky in the CoE plane\nCoE at origin (left), us at distance D; color = l=2 signal around our visible bubble')
ax.grid(alpha=.3)
plt.tight_layout();plt.savefig('quad_xy_coe.png',dpi=120)
print('a2 (dT/T l=2 coeff)',a2,'  l=2 amp uK:',round(T0*a2,0))
print('q at psi=0 (toward CoE):',round(q[0],0),'uK; psi=90:',round(q[180],0),'; psi=180:',round(q[-1],0))
