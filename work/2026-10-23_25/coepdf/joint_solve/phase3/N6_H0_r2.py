# FIXED 3:3x PM: curve = squared-refit H/<H> (migkas_refit_sq.py), old CSV was H-like so sqrt was double-applied.
# NEW phase3: r2/chi2 of the N6 curves vs Migkas points in H0 space at <H0> = mean/low/high from h0_sources.json (same transform as N6_H0_direction.py)
import numpy as np,json
J=json.load(open('/workspace/coepdf/joint_solve/phase3/h0_sources.json'))
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];el,eh=d['A_err_low'],d['A_err_high']
c=np.loadtxt('/workspace/coepdf/joint_solve/phase3/migkas_refit_sq_curve.csv',delimiter=',',skiprows=1);lc,A,Hc=c[:,0],c[:,1],c[:,2]  # squared refit: A=af^2/<af^2>, Hc=af/<af> (H-like, no extra sqrt)
o={}
for name,H in (('mean',J['H0_mean_invvar']),('low',J['H0_low']),('high',J['H0_high'])):
    h=Hc.copy();s=H/h.mean();k=H/np.mean(y**0.5);hd=y**0.5*k;m=np.interp(l,lc,h*s,period=360)
    e=0.5*k*(((y+eh)**0.5-(y-el).clip(1e-9)**0.5))
    o[name]=dict(H0=H,r2=float(1-np.sum((hd-m)**2)/np.sum((hd-hd.mean())**2)),chi2=float(np.sum(((hd-m)/e)**2)),chi2_flat=float(np.sum(((hd-H)/e)**2)),n=len(l))
json.dump(o,open('/workspace/coepdf/joint_solve/phase3/N6_H0_r2.json','w'),indent=1);print(o)
