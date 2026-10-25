import sys,json,warnings,pickle;warnings.filterwarnings('ignore');sys.path.insert(0,'.');sys.path.insert(0,'..');from static import *
UB=10**0.021927885054513485;FB=0.015426145753569655
obs=np.load('../../../three_vel/point/cmb/obs_ring.npy')[:-1];PR=float(np.ptp(obs))
psis=np.round(np.arange(150,180.0001,0.01),4)
res={}
for tag,f in (('nomass',0.0),('mass',FB)):
    m=Model(UB,f);r=sky(m,psis);res[tag]=dict(m=dict(u=m.u,f=f,M=m.M,t0=m.t0,tc=m.tc),r=r)
    st=np.array([x['status'] for x in r]);opz=np.array([x.get('opz',np.nan) for x in r])
    print(tag,'t0',m.t0,'tc',m.tc,'D (green 4.2 Gpc)',4.2/(m.t0-m.tc))
    prev=None
    for x in r:
        if x['status']!=prev: print('   from psi',x['psi'],x['status'],'1+z %.5g num %.5g den %.5g'%(x['opz'],x['num'],x['den']));prev=x['status']
    print('   1+z at 150 (cap edge) %.4f, 180 (centre) %.4f; emitted before our crossing: %s'%(r[0]['opz'],r[-1]['opz'],sorted(set(x.get('before_our_crossing') for x in r))))
    if tag=='nomass': print('   max |1+z - t0/te| =',np.nanmax([abs(x['opz']-m.t0/x['te']) for x in r]))
    g=np.isfinite(opz)&(opz>0);bad=~g
    print('   directions with 1+z<=0 or infinite: psi',(psis[bad].min(),psis[bad].max()) if bad.any() else None)
    # dT/T over cap (sin-weighted), only defined when no divergence
    p=np.radians(psis);w=np.sin(p);w/=w.sum()
    T=1/opz
    if bad.any():
        res[tag]['stats']=None;print('   dT/T not defined (1+z passes through 0 / infinity)')
    else:
        d=T/np.sum(w*T)-1;X=np.c_[np.ones_like(p),np.cos(p)];sw=np.sqrt(w);c=np.linalg.lstsq(X*sw[:,None],d*sw,rcond=None)[0];rr=d-X@c
        s=dict(rms=float(np.sqrt(np.sum(w*d*d))),rng=float(np.ptp(d)),dip=float(abs(c[1])),rrms=float(np.sqrt(np.sum(w*rr*rr))),rrng=float(np.ptp(rr)),d=d,res=rr)
        res[tag]['stats']=s;print('   cap dT/T: RMS %.5f range %.5f (%.1fx Planck) dipole coef %.5f leftover RMS %.5f range %.5f (%.2fx Planck)'%(s['rms'],s['rng'],s['rng']/PR,s['dip'],s['rrms'],s['rrng'],s['rrng']/PR))
    # also separate finite segments for mass
    if bad.any():
        for lo,hi,nm in ((150,psis[bad].min(),'outer segment'),(psis[bad].max(),180.0001,'inner segment')):
            sel=(psis>=lo)&(psis<hi)&g
            print('   %s psi %.2f-%.2f: 1+z min %.4g max %.4g'%(nm,psis[sel].min(),psis[sel].max(),opz[sel].min(),opz[sel].max()))
# comoving at same speed (full sky)
import os;_c=os.getcwd();os.chdir('..');from sweep_dtt import stats;os.chdir(_c)
for t,fn in (('nomass','../zmig_nomass.npz'),('mass','../zmig_mass.npz')):
    os.chdir('..');s=stats(fn.replace('../',''));os.chdir(_c);res['comoving_'+t]=s;print('comoving',t,'RMS %.5f range %.5f (%.1fx Planck) dipole %.5f leftover range %.5f (%.2fx) 1+z %.4f-%.4f D %.3f green_in %s'%(s['rms'],s['rng'],s['rng']/PR,s['dip'],s['rrng'],s['rrng']/PR,s['zmin'],s['zmax'],s['D'],s['green_in']))
    # comoving restricted to the same cap 150-180 for like-for-like
    d=np.array(s['d']);ps=np.arange(0,181.);sel=ps>=150;p=np.radians(ps[sel]);w=np.sin(p);w/=w.sum()
    z=np.load(fn)['o'][:,1][sel];T=1/z;dc=T/np.sum(w*T)-1;print('   comoving same 150-180 cap: range %.5f (%.1fx Planck) 1+z %.4f-%.4f'%(np.ptp(dc),np.ptp(dc)/PR,z.min(),z.max()))
pickle.dump(dict(res=res,PR=PR,psis=psis),open('best.pkl','wb'))
