import numpy as np,healpy as hp
NS=64;R=21.7439;EPS=0.025*R;n=np.array(hp.pix2vec(NS,np.arange(hp.nside2npix(NS)))).T
def kern(v):
    d=v[None,:]-R*n;dd=(np.sum(d*d,1)+EPS**2)**1.5;k=np.sum(n*d,1)/dd-n@v/np.linalg.norm(v)**3;return k-k.mean()
def shell(r,N=20000):
    i=np.arange(N)+0.5;ph=np.arccos(1-2*i/N);th=np.pi*(1+5**0.5)*i
    return r*np.c_[np.cos(th)*np.sin(ph),np.sin(th)*np.sin(ph),np.cos(ph)]
one=kern(np.array([5.,0,0])).std()
for r in (3.,8.,15.):
    m=sum(kern(p) for p in shell(r))/20000
    print(r,'shell-avg kernel rms',m.std(),'single body rms',one,'ratio',m.std()/one)
