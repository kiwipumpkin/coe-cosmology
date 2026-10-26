# Phase3 (NEW) copy of phase2 p6, output to phase3.
# Phase2 B2: TT chi2 of the seed-42 sky with the compensated low-l sky from p5 run c (actual bodies cancelling (1-f) of field l=2).
# Seed-42 spectrum: l<=40 from Dl_sky (anafast of Resid+Fmap map, ns64), l>40 Dfix+Df*g42. Ensemble: Df+Dfix (bodies change one sky, not the ensemble).
# Usage: python3 p6_tt_compensated.py FIELD_NPZ(with Df,Dfix or Dfix_new) P5_NPZ TAG
import sys,os,json,numpy as np
O='/workspace/coepdf/joint_solve/phase4/';fn,p5,tag=sys.argv[1],sys.argv[2],sys.argv[3]
os.chdir('/workspace/redweb2/stochfield');sys.path.insert(0,'.')
import acoustic_fit as A
z=np.load(fn);Df=z['Df'];Dfix=z['Dfix'] if 'Dfix' in z else z['Dfix_new'];c=np.load(p5);Ds=c['Dl_sky']
L=np.arange(len(Df));g=np.random.default_rng(42).chisquare(2*np.maximum(L,1)+1)/(2*np.maximum(L,1)+1)
D=Df+Dfix;Dr=Dfix+Df*g;Drc=Dr.copy();n=41;Drc[2:n]=Ds[2:n]   # only l<=40: ns64 map spectra unreliable above ~2NS; bodies targeted l<=40
rg=dict(l2_30=(2,30),l30_2500=(30,2500),l2_2500=(2,2500))
out=dict(field=fn,p5=p5,check_l3_10_Dlsky_over_Dr=(Ds[3:11]/Dr[3:11]).tolist(),D2=dict(ens=float(D[2]),seed42=float(Dr[2]),seed42_comp=float(Drc[2])),
 mean={k:A.stats(D,Df,*r) for k,r in rg.items()},seed42={k:A.stats(Dr,Df,*r) for k,r in rg.items()},seed42_comp={k:A.stats(Drc,Df,*r) for k,r in rg.items()})
json.dump(out,open(O+f'p6_{tag}.json','w'),indent=1,default=float);np.savez(O+f'p6_{tag}.npz',D=D,Dr=Dr,Drc=Drc,Df=Df,Lt=A.Lt,Pf=A.Pf)
for s in('mean','seed42','seed42_comp'):print(tag,s,{k:round(out[s][k]['chi2_planck'],2) for k in rg},flush=True)
print('D2',out['D2'],'Dl_sky/Dr l3-10',np.round(out['check_l3_10_Dlsky_over_Dr'],3))
