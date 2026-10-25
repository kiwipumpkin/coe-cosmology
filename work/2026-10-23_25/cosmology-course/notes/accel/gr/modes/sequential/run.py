import json,pickle,numpy as np,multiprocessing as mp;from seq import *
TH=np.arange(0,360,2.0)
CASES=[('temp',0.1),('temp',0.05),('temp',0.01),('temp',0.003)]
def job(a):
    p,(rule,par)=a;u=U(p);ev=ev_temp(par) if rule=='temp' else ev_front(0.05*T1,par);out=[]
    for th in np.radians(TH):
        try:out.append(trace(u,th,ev))
        except Exception as e:out.append(None)
    return a,out
if __name__=='__main__':
    jobs=[(p,c) for p in (1.0,1.08) for c in CASES]
    with mp.Pool(8) as P:R=dict(P.map(job,jobs))
    pickle.dump(R,open('seq.pkl','wb'))
    rows=[]
    for (p,(rule,par)),out in R.items():
        ok=np.array([o is not None for o in out]);T=np.array([1/o['opz'] if o else np.nan for o in out])
        d=T/np.nanmean(T)-1;dd=np.nan_to_num(d);m=np.fft.rfft(dd)/len(dd)
        rows.append(dict(p=p,rule=rule,par=par,missing=int((~ok).sum()),rms=float(np.sqrt(np.nanmean(d**2))),range=float(np.nanmax(d)-np.nanmin(d)),
            dipole=float(2*abs(m[1])),dipole_sign=float(np.sign(m[1].real)),quad=float(2*abs(m[2])),
            te_min=float(np.nanmin([o['te'] for o in out if o])),te_max=float(np.nanmax([o['te'] for o in out if o])),
            mean_1pz=float(np.nanmean([o['opz'] for o in out if o]))))
        print(rows[-1])
    json.dump(rows,open('seq_stats.json','w'),indent=1)
