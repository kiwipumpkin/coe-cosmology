import sys,os;sys.argv=['x','exp']
src=open('unified2.py').read().split("if __name__=='__main__'")[0]
# old best-fit: evidence exp_free maxL, p=1; convert mc->M at H70 (t0 from H70)
import time
for ring in ("low",):
    os.environ['RING']=ring; ns={}; exec(src,ns)
    D,mc,H=40.6904,0.01,70.0
    M=mc*((H*ns['HU'])**2*(D*ns['GPCm'])**3)/(2*ns['G'])
    e=ns['model'](dict(D=D,M=M,l0=15.336,L=1.23538,H=H,qc=0.23,p=1.0))
    print(ring,'chi2',round(e['chi2'],2),'R2',round(e['R2'],3),'M %.3e'%M)
