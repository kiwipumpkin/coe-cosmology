import pickle,functools,numpy as np,sys
sys.argv=['x']
import run_exactc as X
src=open('fill_exactc.py').read().split('d=pickle.load')[0]; g={}; exec(src,g)
g['LTB']=functools.partial(X.LTB,Rdot=6.0,rmin=7.9)
d=pickle.load(open('recession6_fixedmass.pkl','rb')); k=(0.0946/4,0.00123); dt,x,t1=d[k]
d[k]=(dt,g['runfill'](k[0],k[1],t1+dt,x),t1); print(np.isnan(d[k][1]).sum()); pickle.dump(d,open('recession6_fixedmass.pkl','wb'))
