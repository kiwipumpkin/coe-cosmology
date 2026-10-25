import json,pickle,numpy as np,multiprocessing as mp;from rings import *
TOBS=3.0;CASES=[(1.6,12.0),(0.8,2.0),(1.6,0.0),(0.8,0.0)];C=[20.0,50.0];TREC=[2.85,2.9,2.95]
def job(a):
    (k,G),c,tr=a;F=Flow(k,G);ar=float(F.aloc(3.0,tr));return a,dict(arec=ar,res=sky(F,c,TOBS,ar,galilean=True))
if __name__=='__main__':
    jobs=[(cs,c,tr) for cs in CASES for c in C for tr in TREC]
    with mp.Pool(8) as P:R=dict(P.map(job,jobs))
    pickle.dump(R,open('rings_gal.pkl','wb'));rows=[]
    for ((k,G),c,tr),d in R.items():
        x=d['res'];T=x[:,0];ok=np.isfinite(T);dd=T/np.nanmean(T)-1;z=np.nan_to_num(dd);m=np.fft.rfft(z)/len(z)
        rows.append(dict(k=k,GM=G,c=c,t_rec_local=tr,a_rec=d['arec'],missing=int((~ok).sum()),rms=float(np.sqrt(np.nanmean(dd**2))),range=float(np.nanmax(dd)-np.nanmin(dd)),
          dipole=float(2*abs(m[1])),dipole_sign=float(np.sign(m[1].real)),quad=float(2*abs(m[2])),te_0=float(x[0,1]),te_90=float(x[45,1]),te_180=float(x[90,1]),mean_T_ratio=float(np.nanmean(T))))
        print({kk:(round(v,6) if isinstance(v,float) else v) for kk,v in rows[-1].items()})
    json.dump(rows,open('rings_stats_gal.json','w'),indent=1)
