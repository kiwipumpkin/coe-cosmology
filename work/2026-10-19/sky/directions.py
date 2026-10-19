import numpy as np,json
from astropy.coordinates import SkyCoord;import astropy.units as u
def vec(l,b): l,b=np.radians(l),np.radians(b);return np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
def sep(a,b): return float(np.degrees(np.arccos(np.clip(vec(*a)@vec(*b),-1,1))))
def axsep(a,b): s=sep(a,b);return min(s,180-s)
L0=15.61   # Migkas l0 (point fits 15.59-15.65; profile 15.605-15.612; soft fixed 15.2)
D=dict(
 CoE=(L0+180,0.0), antiCoE=(L0,0.0),
 Planck_dipole_hot=(264.021,48.253), Planck_dipole_cold=(84.021,-48.253),
 quadrupole_axis_dOC04=(257.73,57.40), octopole_axis_dOC04=(238.50,60.44),
 CatWISE_dipole=(238.2,28.8),
 Migkas21_H0_low=(280.0,-15.0), Migkas21_LxT_bulkflow_MR=(315.0,-10.0), Migkas21_YszT_bulkflow_MA=(267.0,-28.0))
out={}
for k,(l,b) in D.items():
    c=SkyCoord(l=l*u.deg,b=b*u.deg,frame='galactic').icrs;out[k]=dict(l=l,b=b,ra=float(c.ra.deg),dec=float(c.dec.deg))
    print(f"{k:28s} l {l:7.2f} b {b:6.2f}  RA {c.ra.deg:7.2f} Dec {c.dec.deg:6.2f}  sep from CoE {sep((l,b),D['CoE']):6.1f}  from anti-CoE {sep((l,b),D['antiCoE']):6.1f}  axis sep {axsep((l,b),D['CoE']):5.1f}")
# CoE allowed ring: longitude fixed, latitude free -> min separation of each reference from the meridian l=L0+180 (and l=L0)
for k in ('Planck_dipole_hot','quadrupole_axis_dOC04','CatWISE_dipole','Migkas21_H0_low'):
    bs=np.linspace(-90,90,3601);s=[sep(D[k],(L0+180,b)) for b in bs];i=int(np.argmin(s));s2=[sep(D[k],(L0,b)) for b in bs];j=int(np.argmin(s2))
    print(f"{k}: closest point on CoE meridian l={L0+180:.2f}: b={bs[i]:.1f}, sep {s[i]:.1f}; on l={L0:.2f}: b={bs[j]:.1f}, sep {s2[j]:.1f}")
print('Migkas21 low-H0 dir vs CoE axis longitude difference', (280-L0)%360)
json.dump(out,open('directions.json','w'),indent=1)
