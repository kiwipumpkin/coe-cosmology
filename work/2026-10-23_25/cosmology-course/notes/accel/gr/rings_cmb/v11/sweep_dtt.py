import numpy as np,json
obs=np.load('../../three_vel/point/cmb/obs_ring.npy')[:-1];PR=float(np.ptp(obs))
psi=np.arange(0,181,1.0);p=np.radians(psi);w=np.sin(p);w[0]=w[-1]=np.sin(np.radians(0.5))/4;w=w/w.sum()
def stats(fn):
    d=np.load(fn);o=d['o'];z=o[:,1];T=1/z;dd=T/np.sum(w*T)-1
    X=np.c_[np.ones_like(p),np.cos(p)];sw=np.sqrt(w);c=np.linalg.lstsq(X*sw[:,None],dd*sw,rcond=None)[0];res=dd-X@c
    t0,tc,xo=float(d['t0']),float(d['tc']),float(d['xo']);rv=t0-tc;D=4.2/rv
    dt=np.hypot(d['rt']*np.cos(d['phi'])-xo,d['rt']*np.sin(d['phi']))
    return dict(rms=float(np.sqrt(np.sum(w*dd*dd))),rng=float(np.ptp(dd)),dip=float(abs(c[1])),rrng=float(np.ptp(res)),rrms=float(np.sqrt(np.sum(w*res*res))),
      D=D,age=t0*D*3.26156,since=rv*D*3.26156,green_in=bool(np.nanmin(dt)>=rv),gap=float((np.nanmin(dt)-rv)*D),zmin=float(z.min()),zmax=float(z.max()),unres=int((~np.isfinite(z)).sum()),d=dd.tolist())
if __name__=='__main__':
    print('check zc25:',{k:v for k,v in stats('zc25_nomass.npz').items() if k!='d'},'\n',{k:v for k,v in stats('zc25_mass.npz').items() if k!='d'})
    P=json.load(open('../../three_vel/point/prof_c4.json'));out={}
    for k in P:
        a=stats(f'sweep/s{k}_nomass.npz');b=stats(f'sweep/s{k}_mass.npz');out[k]=dict(frac=P[k]['x'][0],R2=P[k]['R2'],chi2=P[k]['chi2'],nomass=a,mass=b,dRMS=100*(b['rms']/a['rms']-1),dRng=100*(b['rng']/a['rng']-1))
    json.dump(dict(PR=PR,rows=out),open('sweep/sweep_dtt.json','w'))
    print('Planck',PR)
    for k,v in out.items():
        for t in ('nomass','mass'):
            s=v[t];print(f"{k:>3}c {t:6} D {s['D']:.3f} age {s['age']:.3f} since {s['since']:.3f} RMS {s['rms']:.5f} range {s['rng']:.5f} r/P {s['rng']/PR:.2f} dip {s['dip']:.5f} resid r/P {s['rrng']/PR:.3f} green_in {s['green_in']} gap {s['gap']:.3f} 1+z {s['zmin']:.4f}-{s['zmax']:.4f}")
        print(f"     mass frac {v['frac']:.6f}  dRMS {v['dRMS']:+.2f}%  dRange {v['dRng']:+.2f}%")
