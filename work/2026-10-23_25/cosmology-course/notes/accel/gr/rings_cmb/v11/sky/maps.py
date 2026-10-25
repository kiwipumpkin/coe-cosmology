import numpy as np,json,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from directions import D,vec,sep
import os
B=os.environ.get('COEB')
LB=os.environ.get('COE')
if LB:
    cl,cb=map(float,LB.split(','));D['CoE']=(cl,cb);D['antiCoE']=((cl+180)%360,-cb)
if B:
    B=float(B);D['CoE']=(195.61,B);D['antiCoE']=(15.61,-B)
if os.environ.get('ROT'):
    D['CoE']=(240.0,30.0);D['antiCoE']=(60.0,-30.0)
obs=np.load('../../../three_vel/point/cmb/obs_ring.npy')[:-1];PR=float(np.ptp(obs))
psi=np.arange(0,181,1.0);p=np.radians(psi);w=np.sin(p);w[0]=w[-1]=np.sin(np.radians(0.5))/4;w/=w.sum()
def prof(fn):
    z=np.load(fn)['o'][:,1];T=1/z;d=T/np.sum(w*T)-1;return d
RUNS=[('../distributed/z4_dist_a_own.npz','own crossing, distributed mass, normalisation a (equal mass)'),
      ('../distributed/z4_dist_b_own.npz','own crossing, distributed mass, normalisation b (Migkas-fitted)'),
      ('../distributed/z4_dist_b_sim.npz','simultaneous, distributed mass, normalisation b (Migkas-fitted)'),
      ('../z4_m02.npz','own crossing, point mass, fraction 0.2 (reference)')]
nside=64;pix=np.arange(hp.nside2npix(nside));th,ph=hp.pix2ang(nside,pix);V=np.c_[np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)]
aC=vec(*D['antiCoE']);ang=np.degrees(np.arccos(np.clip(V@aC,-1,1)))
fig=plt.figure(figsize=(18,11.5));res={}
mk=[('CoE','k','*',420,(f'CoE (l = {cl:g}°, b = {cb:+g}°, paper main.tex l.359–362)' if LB else f'CoE (l = 195.6°, b = {B:+.1f}° assumed: best mean match)' if B else 'CoE (l = 240°, b = +30°, Chris rotation)') if (os.environ.get('ROT') or B or LB) else 'CoE (l = 195.6°, b = 0° assumed)'),('antiCoE','k','X',160,(f'anti-CoE (l = {(cl+180)%360:g}°, b = {-cb:+g}°)' if LB else f'anti-CoE (l = 15.6°, b = {-B:+.1f}°)' if B else 'anti-CoE (l = 60°, b = −30°)') if (os.environ.get('ROT') or B or LB) else 'anti-CoE (l = 15.6°, b = 0°)'),('Planck_dipole_hot','r','^',140,'Planck dipole hot (264.0°, +48.3°)'),('Planck_dipole_cold','b','v',140,'Planck dipole cold (84.0°, −48.3°)'),
    ('quadrupole_axis_dOC04','m','D',110,'quadrupole axis, de Oliveira-Costa+04 (257.7°, +57.4°)'),('CatWISE_dipole','g','s',110,'CatWISE quasar dipole (238.2°, +28.8°)'),('Migkas21_H0_low','orange','o',120,'Migkas+21 low-H₀ direction (280°, −15°)')]
for i,(fn,lab) in enumerate(RUNS):
    d=prof(fn);m=np.interp(ang,psi,d)
    imax,imin=int(np.argmax(d)),int(np.argmin(d))
    # quadrupole axis of the map (axisymmetric -> along CoE axis); verify numerically
    a2=hp.map2alm(m,lmax=2);a2[hp.Alm.getidx(2,0,0)]=0;a2[hp.Alm.getidx(2,1,0)]=0;a2[hp.Alm.getidx(2,1,1)]=0
    q=hp.alm2map(a2,nside);j=int(np.argmax(np.abs(q)));qa=np.degrees(hp.pix2ang(nside,j,lonlat=True))
    qa=hp.pix2ang(nside,j,lonlat=True)
    res[lab]=dict(range_over_Planck=float(np.ptp(d)/PR),hot_angle_from_antiCoE=float(psi[imax]),cold_angle_from_antiCoE=float(psi[imin]),quad_axis_lb=[float(qa[0]),float(qa[1])],quad_sign_at_axis='hot' if q[j]>0 else 'cold')
    print(lab,res[lab])
    hp.mollview(m,sub=(2,2,i+1),fig=fig.number,title=f"{lab}\nrange {np.ptp(d)/PR:.1f}× Planck ring; hottest {psi[imax]:.0f}° from anti-CoE, coldest {psi[imin]:.0f}°",cmap='RdBu_r',unit='dT/T',min=-np.abs(d).max(),max=np.abs(d).max(),notext=False)
    hp.graticule(dpar=30,dmer=60,alpha=.3)
    bb=np.linspace(-89,89,100)
    for lm in (() if (os.environ.get('ROT') or LB) else (195.61,15.61)): hp.projplot(np.full_like(bb,lm),bb,'k--',lonlat=True,lw=1)
    for k,c,s,sz,_ in mk: hp.projscatter(D[k][0],D[k][1],lonlat=True,color=c,marker=s,s=sz,edgecolors='w',linewidths=.8,zorder=6 if k=='CoE' else 5)
h=[plt.Line2D([],[],color=c,marker=s,ls='',ms=11,markeredgecolor='w') for _,c,s,_,_ in mk]+([] if (os.environ.get('ROT') or LB) else [plt.Line2D([],[],color='k',ls='--')])
fig.legend(h,[t for *_,t in mk]+([] if (os.environ.get('ROT') or LB) else ['dashed: meridians l = 195.6° and 15.6° (CoE latitude not fixed by the Migkas longitude fit)']),loc='upper center',ncol=3,fontsize=11,frameon=False,bbox_to_anchor=(0.5,1.07))
fig.text(0.5,1.085,'Model dT/T skies in galactic coordinates (Mollweide, l = 0 at centre, l increasing to the left); speed 4c; '+((f'CoE from paper (ST CoE paper 2026-09-25, main.tex lines 359–362): ({cl:g}°, {cb:+g}°)') if LB else (f'CoE l = 195.6° from Migkas ℓ₀+180°; b = {B:+.1f}° is an assumption (minimises mean separation from Planck hot, quadrupole axis, CatWISE)') if B else 'CoE axis from Chris rotation: +x → (l, b) = (−120°, +30°) = (240°, +30°); CoE star sits 2° from CatWISE square' if os.environ.get('ROT') else 'CoE placed at Migkas axis longitude, latitude 0°'),ha='center',fontsize=12)
fig.subplots_adjust(top=0.88,hspace=0.25)
fig.savefig(('../sky_dTT_mollweide_4c_paper.png' if LB else '../sky_dTT_mollweide_4c_l180.png' if B else '../sky_dTT_mollweide_4c_rot.png' if os.environ.get('ROT') else '../sky_dTT_mollweide_4c.png'),dpi=90,bbox_inches='tight')
json.dump(res,open(('model_sky_paper.json' if LB else 'model_sky_l180.json' if B else 'model_sky_rot.json' if os.environ.get('ROT') else 'model_sky.json'),'w'),indent=1)
