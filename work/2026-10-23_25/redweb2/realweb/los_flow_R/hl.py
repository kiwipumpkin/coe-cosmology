import sys;sys.argv=['x'];src=open('flowR.py').read().split('\nif __name__')[0];exec(src)
k=113;Rf,Rq=grids(500000+k);Rq=Rq+tocoarse(Rf)
fq=phi(Gq,Rq,Gq['dx']);ff_=phi(Gf,Rf,Gf['dx']);fs=phi(Gf,Rf,Gq['dx']);ff=[a-b for a,b in zip(ff_,fs)];del ff_,fs
sm32=np.load('/workspace/redweb2/sm_red.npy');asm=hp.map2alm(sm32,lmax=95)
for ns in (32,64,128,256):
    t=time.time();V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T;o=comp(fq,ff,4.0);o['sm']=hp.alm2map(asm,ns,lmax=95)
    np.savez(O+f'hl_R4_ns{ns}.npz',**o);print(ns,'%.0fs'%(time.time()-t),flush=True)
