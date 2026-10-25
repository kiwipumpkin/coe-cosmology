import os,sys,time,json,numpy as np
law,mode,vm=sys.argv[1:4]; os.environ['R0DENSE']='1'
if vm=='gtc': os.environ['VMODE']='gtc'
sys.argv=['x',law,mode]; exec(open('fit_lc.py').read().split("ck=f'fwd")[0]); import grltb as G
x=np.array(json.load(open(f'gr/ck_{law}_{mode}_{vm}.json'))['x']); t=time.time()
P,o=build(x); f=G.fields(o); tr=G.trace(f,np.arange(0,360,4.),nstep=320,lbmax_frac=0.999); af=G.migkas(tr); print(F.chi2_l0(af),time.time()-t, nm, BX if mode=='back' else [B[n] for n in nm], F.CHI0)
