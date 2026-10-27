# Phase5 (NEW): CMB queue (TAG D mc eta L) from Migkas fits: A_d<d> = fitd_<d>.json (best at fixed d, k free -> runs away),
# K<lk> = phase4b profile_k4b A row (d free, k fixed), J_<d>_<lk> = fitdk_<d>_<lk>.json (d,k fixed). Priority order.
import json,glob,os,sys
O='/workspace/coepdf/joint_solve/phase5/';sys.path.insert(0,'/workspace/coepdf/joint_solve/phase4b')
P=json.load(open('/workspace/coepdf/joint_solve/phase4b/profile_k4b.json'));Q=[]
def add(t,r):Q.append(f"{t} {r['D']:.6f} {r['mc']:.6f} {r['eta']:.6f} {r['L']:.6f}")
for d in ['0.001','2','10']:add('A_d'+d,json.load(open(O+f'fitd_{d}.json')))
KA={round(r[1],2):r[3] for r in P if r[0]=='A'}
for lk in (-1.0,1.0,3.0):add('K%g'%lk,KA[lk])
for f in sorted(glob.glob(O+'fitdk_*[0-9].json')):r=json.load(open(f));add('J_%g_%g'%(r['d'],round(__import__('math').log10(r['k']),3)),r)
for d in ['0.5','1','3','4','6','8']:add('A_d'+d,json.load(open(O+f'fitd_{d}.json')))
for lk in (0.0,2.0,5.0):add('K%g'%lk,KA[lk])
print('\n'.join(Q))
