"""Periodic (Ewald) lattice energies, HCP vs FCC, with our sign convention U = -G sum m_i m_j / r_ij (voids = negative mass)
plus a uniform background that cancels the net mass (the 'compensation picture': the mean density is subtracted, as in the
CMB/peculiar-potential calculation). u = U per node / (G M_cluster^2 / d), d = nearest-neighbour spacing 120.9 Mpc.
Gravity is Coulomb with the sign flipped (like masses attract), so u_grav = -E_Coulomb(q = m/M_cluster).
Check: unit charges on FCC/HCP must reproduce the Madelung constants 1.791747 / 1.791676 (per particle, units e^2/a_ws)."""
import json,numpy as np
from scipy.special import erfc
rng=np.random.default_rng(11)
def ewald(X,q,Lbox,alpha=None,kmax=None):
    N=len(X); V=np.prod(Lbox); alpha=alpha or 5.6/Lbox.min(); 
    # real space, 27 images
    Er=0.0; sh=np.array([[i,j,k] for i in (-1,0,1) for j in (-1,0,1) for k in (-1,0,1)])*Lbox
    for s in sh:
        D=np.linalg.norm(X[:,None,:]-X[None,:,:]+s,axis=2)
        if not s.any(): np.fill_diagonal(D,np.inf)
        Er+=0.5*np.sum(np.outer(q,q)*erfc(alpha*D)/D)
    # reciprocal
    nk=np.ceil(kmax or 3.2*alpha*Lbox.max()/np.pi).astype(int)
    g=np.array([[i,j,k] for i in range(-nk,nk+1) for j in range(-nk,nk+1) for k in range(-nk,nk+1) if (i,j,k)!=(0,0,0)])
    K=2*np.pi*g/Lbox; k2=(K**2).sum(1); m=k2<(2*alpha*3.2)**2; K,k2=K[m],k2[m]
    Ek=0.0
    for c in range(0,len(K),4000):
        S=(q[None,:]*np.exp(1j*K[c:c+4000]@X.T)).sum(1); Ek+=np.sum(np.exp(-k2[c:c+4000]/(4*alpha**2))/k2[c:c+4000]*np.abs(S)**2)
    Ek*=2*np.pi/V
    Es=-alpha/np.sqrt(np.pi)*np.sum(q**2); Eb=-np.pi*q.sum()**2/(2*V*alpha**2)
    return Er+Ek+Es+Eb
d=1.0
def fcc_cell(n):
    a=np.sqrt(2)*d; B=np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]]); T=np.array([[i,j,k] for i in (.25,.75) for j in (.25,.75) for k in (.25,.75)])
    O=np.array([[.5,.5,.5],[.5,0,0],[0,.5,0],[0,0,.5]]); g=np.array([[i,j,k] for i in range(n) for j in range(n) for k in range(n)],float)
    f=lambda U:(g[:,None,:]+U[None]).reshape(-1,3)*a; return f(B),f(T),f(O),np.full(3,n*a)
def hcp_cell(n1,n2,n3):
    """orthohexagonal cell (d, sqrt3 d, c), 4 sites, 8 tetra holes, 4 octa holes"""
    c=np.sqrt(8/3)*d; Lc=np.array([d,np.sqrt(3)*d,c])
    A2=np.array([[0,0,0],[.5,.5,0]]); Bs=A2+[.5,1/6,.5]    # A layer and B layer (B over the 'up' triangles)
    S=np.vstack([A2,Bs%1]); 
    # holes: tetra above A sites at 3/8 c and below at -3/8 c (=5/8), above/below B sites; octa at the C positions
    Tt=np.vstack([A2+[0,0,3/8],A2+[0,0,5/8],Bs+[0,0,-3/8],Bs+[0,0,3/8]])%1; Cp=A2+[.5,5/6,0]; Oo=np.vstack([Cp+[0,0,.25],Cp+[0,0,.75]])%1
    g=np.array([[i,j,k] for i in range(n1) for j in range(n2) for k in range(n3)],float)
    f=lambda U:((g[:,None,:]+U[None])*Lc).reshape(-1,3); return f(S),f(Tt),f(Oo),Lc*[n1,n2,n3]
def check_holes(C,T,O,L):
    from scipy.spatial import cKDTree
    t=cKDTree(C%L,boxsize=L); a=t.query(T%L,5)[0]; b=t.query(O%L,7)[0]
    return np.round(a[:,:4].max(),4),np.round(a[:,4].min(),4),np.round(b[:,:6].max(),4),np.round(b[:,6].min(),4)
FC=fcc_cell(4); HC=hcp_cell(6,4,4)   # 256 and 384 clusters, boxes ~5.7 d and 6.0 x 6.9 x 6.5 d
out={"check_holes (tetra 4th, 5th / octa 6th, 7th nearest cluster, units d)":{"FCC":check_holes(*FC),"HCP":check_holes(*HC)}}
print(out)
aws=lambda N,V:(3*V/(4*np.pi*N))**(1/3)
for nm,(C,T,O,L) in (("FCC",FC),("HCP",HC)):
    E=ewald(C,np.ones(len(C)),L); out.setdefault("Madelung check (unit charges, -E/N * a_ws)",{})[nm]=float(-E/len(C)*aws(len(C),np.prod(L)))
print(out)
MC,MV=1.40e50,-8.39e50; mv=MV/MC
def u(C,H,L,fill,qv,nrep=1):
    vals=[]
    for r in range(nrep):
        Hs=H[rng.random(len(H))<fill] if fill<1 else H
        X=np.vstack([C,Hs]); q=np.r_[np.ones(len(C)),np.full(len(Hs),qv)]
        vals.append(-ewald(X,q,L)/len(X))
    return float(np.mean(vals)),float(np.std(vals)/np.sqrt(max(nrep-1,1))) if nrep>1 else 0.0
R={}
cases=[("clusters only",None,0,0,1),
       ("web ratio: 1/12 of tetra holes, real masses (void = -5.99 cluster)","T",1/12,mv,40),
       ("web ratio: 1/6 of octa holes, real masses","O",1/6,mv,40),
       ("all tetra holes, compensated (void = -0.5 cluster)","T",1,-0.5,1),
       ("all octa holes, compensated (void = -1 cluster; FCC = rock-salt, HCP = NiAs pattern)","O",1,-1.0,1),
       ("all tetra holes, real masses + background","T",1,mv,1),
       ("all octa holes, real masses + background","O",1,mv,1)]
for nm,h,fill,qv,nrep in cases:
    R[nm]={}
    for lat,(C,T,O,L) in (("FCC",FC),("HCP",HC)):
        if h is None: R[nm][lat]=(float(-ewald(C,np.ones(len(C)),L)/len(C)),0.0)
        else: R[nm][lat]=u(C,T if h=="T" else O,L,fill,qv,nrep)
    f,hh=R[nm]["FCC"][0],R[nm]["HCP"][0]; R[nm]["lower energy (more bound)"]="HCP" if hh<f else "FCC"; R[nm]["HCP - FCC"]=hh-f
    print(nm,R[nm],flush=True)
out["u per node, units G M_cl^2/d, with background (mean ± s.e. over random fills)"]=R
json.dump(out,open("hcp_lattice/hcp_ewald.json","w"),indent=1)
