import sys;sys.argv=['x']
src=open('/workspace/redweb2/realweb/los_flow_R/flowR.py').read().split('\nif __name__')[0];exec(src)
k=113;Rf,Rq=grids(500000+k);Rq=Rq+tocoarse(Rf)
fq=phi(Gq,Rq,Gq['dx']);ff_=phi(Gf,Rf,Gf['dx']);fs=phi(Gf,Rf,Gq['dx']);ff=[a-b for a,b in zip(ff_,fs)];del ff_,fs
res={}
for R in np.arange(2,20.01,0.5):
    t=time.time();o=comp(fq,ff,R)
    for q,v in o.items():res[f'R{R:g}_{q}']=v
    print(R,'%.1fs'%(time.time()-t),flush=True)
    np.savez(O+'u113_Rscan.npz',**res)
