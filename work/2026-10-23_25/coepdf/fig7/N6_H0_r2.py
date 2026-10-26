# r2 and chi2 of the red <H0>=70 curve vs Migkas points in H0 space (same transform as N6_H0_direction.py)
import numpy as np,json
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];el,eh=d['A_err_low'],d['A_err_high']
c=np.loadtxt('/workspace/cosmology-course/notes/catwise-eta-migkas-fit-d14.3.csv',delimiter=',',skiprows=1);lc,A=c[:,0],c[:,1]
H=70.;h=A**0.5;s=H/h.mean();k=H/np.mean(y**0.5)
hd=y**0.5*k;m=np.interp(l,lc,h*s,period=360)
e=0.5*k*(((y+eh)**0.5-(y-el).clip(1e-9)**0.5))
r2=1-np.sum((hd-m)**2)/np.sum((hd-hd.mean())**2);chi2=np.sum(((hd-m)/e)**2);chi2n=np.sum(((hd-H)/e)**2)
o=dict(r2=r2,chi2=chi2,chi2_flat70=chi2n,n=len(l));json.dump(o,open('N6_H0_r2.json','w'),indent=1);print(o)
