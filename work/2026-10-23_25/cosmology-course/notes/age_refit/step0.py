import os,json,numpy as np
src=open('unified2.py').read().split("if __name__=='__main__'")[0]; out={}
import sys; sys.argv=['x','exp']
for ring in ('old','low'):
    os.environ['RING']=ring; ns={}; exec(src,ns)
    D,mc,H=40.6904,0.01,70.0; M=mc*((H*ns['HU'])**2*(D*ns['GPCm'])**3)/(2*ns['G'])
    e=ns['model'](dict(D=D,M=M,l0=15.336,L=1.23538,H=H,qc=0.3,p=1.0),True); out[ring]=dict(chi2=e['chi2'],R2=e['R2'],A=e['A'])
json.dump(out,open('step0.json','w'));print({k:(v['chi2'],v['R2']) for k,v in out.items()})
