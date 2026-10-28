from astropy.io import fits; from astropy.coordinates import SkyCoord; import astropy.units as u, numpy as np, pandas as pd
for f,z in [('erass1cl_primary_v3.2.fits','BEST_Z'),('erass1cl_cosmology_v1.1.fits','Z_LAMBDA')]:
    d=fits.open(f)[1].data; c=SkyCoord(d['RA']*u.deg,d['DEC']*u.deg).galactic
    cols={'NAME':d['NAME'],'RA':d['RA'],'DEC':d['DEC'],'GLON':c.l.deg,'GLAT':c.b.deg,'z':d[z]}
    for k in ['L500','L500_L','L500_H','KT','KT_L','KT_H','YX500','YX500_L','YX500_H']: cols[k]=d[k]
    pd.DataFrame({k:np.asarray(v).astype(str) if k=='NAME' else np.asarray(v,dtype=float) for k,v in cols.items()}).to_csv(f.replace('.fits','_lbz.csv'),index=False)
    print(f, 'KT finite:', np.isfinite(d['KT']).sum(), 'KT>0:', (d['KT']>0).sum())
