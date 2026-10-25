import pickle,json,glob,numpy as np,sys; sys.path.insert(0,'..'); from grcore import stats,norm
M=pickle.load(open('mig.pkl','rb'))
for f in sorted(glob.glob('cmb_*.pkl')):
    d=pickle.load(open(f,'rb')); x=d['x']; nn=int(np.isnan(x).sum())
    if nn: x=x.copy(); b=np.isnan(x); x[b]=np.interp(np.where(b)[0],np.where(~b)[0],x[~b])
    x=norm(x); s=stats(x); m=M.get((d['i'],d['C']),{})
    print(d['i'],'%+.1f'%d['C'],'missing',nn,'dip %.5f quad %.5f cold %.5f rms %.5f'%(s[0],s[1],x.min(),np.sqrt(np.mean(x[:-1]**2))),'mig chi2 %.3f R2 %.4f'%(m.get('chi2',np.nan),m.get('R2',np.nan)))
