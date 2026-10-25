import sys,json,numpy as np;sys.argv=['x','mass','m345','0']
src=open('lvn345.py').read().split('fam,which,NM=')[0];exec(src)
out={'start':json.load(open('start.json'))}
c=cs;out['start_1e5']={f'l{l}':float(lrms(c,[l])/1e-5) for l in range(2,13)};out['start_1e5']['l6_12']=float(lrms(c,list(range(6,13)))/1e-5)
for f in ('mass','void','both'):
    R=json.load(open(f'res_{f}_m345.json'));c=total(R['params'],R['types'])
    T=TYP;sz=[];
    for p,t in zip(R['params'],R['types']):sz.append(dict(type=t,cap_deg=float(np.degrees(p[2]*T[t]['rmax']/DLSS)),contrast=float(T[t]['dl']+p[3]*(T[t]['dh']-T[t]['dl']))))
    out[f]=dict(curve=R['curve'],first_N_1x=R['first_N_1x'],best=float(min(v for n,v in R['curve'])),bodies=sz,
      after_1e5={'l2':float(lrms(c,[2])/1e-5),'l345':float(lrms(c,[3,4,5])/1e-5),'l6_12':float(lrms(c,list(range(6,13)))/1e-5)},dipole_after_xPlanck=float(met(c)[0]),v_kms_after=float(met(c)[0]*PD*299792.458))
    print(f,out[f]['first_N_1x'],out[f]['best'],out[f]['after_1e5'],out[f]['dipole_after_xPlanck'])
print(out['start_1e5']);json.dump(out,open('m345_summary.json','w'),indent=1)
