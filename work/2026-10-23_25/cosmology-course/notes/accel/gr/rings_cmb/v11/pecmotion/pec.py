import os
PL=os.environ.get('PLANCK');GL=os.environ.get('GL')
import numpy as np,healpy as hp,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
PR=float(np.ptp(np.load('../../../three_vel/point/cmb/obs_ring.npy')[:-1]))
beta=370/299792.458;gam=1/np.sqrt(1-beta**2)
psi=np.arange(0,181.0);p=np.radians(psi);w=np.sin(p);w[0]=w[-1]=np.sin(np.radians(.5))/4;w/=w.sum()
RUNS=[('../z4_nomass.npz','own crossing, no mass','C0'),('../distributed/z4_dist_b_own.npz','own crossing, distributed, normalisation b (Migkas-fitted)','C1'),
      ('../z4_sim.npz','simultaneous, no mass','C2'),('../distributed/z4_dist_b_sim.npz','simultaneous, distributed, normalisation b','C3')]
DIRS=[0,90,180,270]  # velocity angle from anti-CoE, in plotted plane
nside=128;V=np.array(hp.pix2vec(nside,np.arange(hp.nside2npix(nside)))).T  # anti-CoE = +x, plane = x-y
chi=np.arange(0,360.0,0.5);c=np.radians(chi);Nc=np.c_[np.cos(c),np.sin(c),0*c]
def Tmodel(T1,n):  # n unit vectors in model rest frame; axisymmetric about +x
    return T1(np.degrees(np.arccos(np.clip(n[:,0],-1,1))))
def obs(T1,vdeg,n):
    if vdeg is None: return Tmodel(T1,n)
    a=np.radians(vdeg);v=np.array([np.cos(a),np.sin(a),0.]);mu=n@v
    # rest-frame direction n' (Lorentz aberration, photon arriving from n): n' = [n + ((γ-1)mu - γβ) v]/(γ(1-β mu))
    npr=(n+((gam-1)*mu-gam*beta)[:,None]*v)/(gam*(1-beta*mu))[:,None]
    npr/=np.linalg.norm(npr,axis=1)[:,None]
    return Tmodel(T1,npr)/(gam*(1-beta*mu))

def planck(ax):
    if GL:
        y1=ax.get_ylim()[1]
        for x,col,t in ((264.02,'r','Planck hot l = 264.02°'),(84.02,'b','Planck cold l = 84.02°'),(180,'k','CoE l = 180°'),(0.3,'k','anti-CoE l = 0°'),(359.7,'k','anti-CoE')):
            ax.axvline(x,color=col,lw=1.6 if col!='k' else 1.0,ls='-' if col!='k' else '--',zorder=10);ax.text(x,y1,' '+t,rotation=90,va='top',ha='left' if x<300 else 'right',fontsize=8,color=col)
        return
    if not PL: return
    y1=ax.get_ylim()[1]
    for x,c,t in ((82.8,'r','Planck hot'),(277.2,'r','Planck hot'),(97.2,'b','Planck cold'),(262.8,'b','Planck cold')):
        ax.axvline(x,color=c,lw=1.4);ax.text(x,y1,f'{t} {x}°',rotation=90,va='top',ha='right' if c=='r' else 'left',fontsize=8,color=c)
CAP='Red/blue lines: observed Planck dipole hot/cold ends as angular distances from the anti-CoE (82.8°, 97.2°; mirrored 277.2°, 262.8°) for the paper CoE (180°, −15°); exact positions along the sideways great circles would differ.'
CAPG='Paper CoE (180°, −15°). Red/blue: Planck dipole hot l = 264.02° / cold l = 84.02°; their latitudes (+48.25° / −48.25°) are off this circle, so the lines mark longitude only.\nMotion directions: 0° = toward anti-CoE (0°, +15°), 90° = toward (90°, 0°), 180° = toward CoE (180°, −15°), 270° = toward (270°, 0°).'
def tag(f):
    if GL: return f.replace('.png','_gl.png')
    return f.replace('.png','_planck.png') if PL else f
LG=np.degrees(np.arctan2(np.sin(c),np.cos(c)*np.cos(np.radians(15))))%360
ORD=np.argsort(LG);XS=LG[ORD] if GL else chi
sel=(lambda y: y[ORD]) if GL else (lambda y: y)
res={};curves={}
for fn,name,col in RUNS:
    z=np.load(fn)['o'][:,1];T1=CubicSpline(psi,1/z)
    d0=1/z/np.sum(w/z)-1;X=np.c_[np.ones_like(p),np.cos(p)];sw=np.sqrt(w)
    c0,c1=np.linalg.lstsq(X*sw[:,None],d0*sw,rcond=None)[0]   # fixed axial fit on no-motion sky
    for vd in [None]+DIRS:
        Tm=obs(T1,vd,V).mean()        # full-sky mean of the observed (or model) map
        d=obs(T1,vd,Nc)/Tm-1;r=d-c0-c1*np.cos(c)
        curves[(name,vd)]=(d,r);res.setdefault(name,{})[str(vd)]=dict(range_over_Planck=float(np.ptp(d)/PR),leftover_over_Planck=float(np.ptp(r)/PR))
for n in res:
    b=res[n]['None']
    for k in res[n]:
        res[n][k]['range_change_over_Planck']=res[n][k]['range_over_Planck']-b['range_over_Planck'];res[n][k]['leftover_change_over_Planck']=res[n][k]['leftover_over_Planck']-b['leftover_over_Planck']
    print(n);[print(f"  v at {k:>4}: range {v['range_over_Planck']:.3f}x ({v['range_change_over_Planck']:+.3f}), leftover {v['leftover_over_Planck']:.3f}x ({v['leftover_change_over_Planck']:+.3f})") for k,v in res[n].items()]
json.dump(res,open('pec_stats.json','w'),indent=1)
DL={0:'toward anti-CoE (0°)',90:'sideways (90°)',180:'toward CoE (180°)',270:'sideways (270°)'}
if GL: DL={0:'toward anti-CoE (l, b) = (0°, +15°)',90:'sideways toward (90°, 0°)',180:'toward CoE (180°, −15°)',270:'sideways toward (270°, 0°)'}
for key,j,fname,yl in (('d',0,'pec_dTT_2x2.png','dT/T'),('r',1,'pec_dTT_2x2_residual.png','dT/T minus fixed no-motion axial monopole + cos dipole')):
    fig,axs=plt.subplots(2,2,figsize=(20,13),sharex=True);
    for ax,vd in zip(axs.flat,DIRS):
        for fn,name,col in RUNS:
            s=res[name][str(vd)];y=curves[(name,vd)][j]
            if key=='d': lab=f"{name}: range {s['range_over_Planck']:.3f}× Planck (change {s['range_change_over_Planck']:+.3f})"
            else: lab=f"{name}: leftover {s['leftover_over_Planck']:.3f}× Planck (change {s['leftover_change_over_Planck']:+.3f})"
            ax.plot(XS,sel(y),color=col,lw=2,label=lab)
        ax.axvline(({0:0.3,90:90,180:180,270:270}[vd] if GL else vd),color='k',ls=':',lw=1);planck(ax);ax.set_xlim(0,360);ax.set_xticks(range(0,361,30));ax.grid(alpha=.25);ax.set_ylabel(yl)
        ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),fontsize=10,frameon=False,title=f'observer moving 370 km/s {DL[vd]} (dotted line = motion direction)')
    for ax in axs[1]: ax.set_xlabel('galactic longitude l (°) on great circle through anti-CoE (0°, +15°) and CoE (180°, −15°),\ntilted 15°; passes (90°, 0°), (270°, 0°)' if GL else 'angle along great circle through CoE axis and velocity, from anti-CoE (°); 180° = CoE')
    fig.suptitle(('dT/T with exact aberration + Doppler, speed 4c' if key=='d' else 'Residual: only the no-motion axial dipole (fit on the no-motion sky, held fixed) is removed, so the 370 km/s motion dipole stays visible')+f'; Planck ring range {PR:.7f}',y=1.0,fontsize=13)
    fig.tight_layout()
    if PL or GL: fig.text(0.5,-0.045,CAPG if GL else CAP,ha='center',fontsize=11)
    fig.savefig(tag(fname),dpi=80,bbox_inches='tight');plt.close(fig)
# summary: sim distributed b
name=RUNS[3][1];fig,ax=plt.subplots(figsize=(16,9))
for vd,col in zip(DIRS,['C0','C1','C2','C3']):
    s=res[name][str(vd)];ax.plot(XS,sel(curves[(name,vd)][1]),color=col,lw=2,label=f"370 km/s {DL[vd]}: range {s['range_over_Planck']:.3f}× Planck ({s['range_change_over_Planck']:+.3f}), leftover {s['leftover_over_Planck']:.3f}× ({s['leftover_change_over_Planck']:+.3f})")
b=res[name]['None'];ax.set_xlim(0,360);ax.set_xticks(range(0,361,30));ax.grid(alpha=.25);planck(ax)
if PL or GL: fig.text(0.5,0.0,CAPG if GL else CAP,ha='center',fontsize=10)
ax.set_xlabel('galactic longitude l (°) on great circle through anti-CoE (0°, +15°) and CoE (180°, −15°),\ntilted 15°; passes (90°, 0°), (270°, 0°)' if GL else 'angle along great circle through CoE axis and velocity, from anti-CoE (°); 180° = CoE');ax.set_ylabel('dT/T minus fixed no-motion axial monopole + cos dipole')
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),fontsize=11,frameon=False,title=f"simultaneous, distributed, normalisation b — no motion: range {b['range_over_Planck']:.3f}× Planck, leftover {b['leftover_over_Planck']:.3f}× (curves = 4 motion directions)")
fig.savefig(tag('pec_sim_b_summary.png'),dpi=90,bbox_inches='tight')
