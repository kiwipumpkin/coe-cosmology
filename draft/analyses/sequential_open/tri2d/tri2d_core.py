"""2D triangular crystal, equal numbers of clusters and voids, U = -G sum M_i M_j / r_ij (voids negative mass). Chris, 10 Oct.
Units: lengths in a (lattice spacing), masses in M_c (cluster), energies in G M_c^2 / a."""
import numpy as np, numba
from scipy.special import erfc
S3=np.sqrt(3)/2
def lattice(Lx,Ly):
    assert Ly%2==0
    I,J=np.meshgrid(np.arange(Lx),np.arange(Ly),indexing="xy"); I=I.ravel(); J=J.ravel()
    X=np.c_[I+0.5*(J%2),J*S3]; idx=lambda i,j: (j%Ly)*Lx+(i%Lx)
    nb=[]
    for i,j in zip(I,J):
        o=[(1,0),(-1,0)]+([(-1,1),(0,1),(-1,-1),(0,-1)] if j%2==0 else [(0,1),(1,1),(0,-1),(1,-1)])
        nb.append([idx(i+di,j+dj) for di,dj in o])
    NB=np.array(nb); E=np.unique(np.sort(np.c_[np.repeat(np.arange(len(I)),6),NB.ravel()],1),axis=0)
    return X,I,J,NB,E,(Lx,Ly*S3)
def ewald_phi(D,W,alpha=None,nimg=2,kidx=None):
    """2D-periodic Ewald sum of 1/r (3D Coulomb kernel, sources in the plane) with neutralizing background (k=0 term).
    D: (n,2) displacements. Returns phi(D); for D=0 returns the self-image sum (n=0 term excluded)."""
    Wx,Wy=W; A=Wx*Wy; alpha=alpha or 6.0/min(W); out=np.zeros(len(D))
    for nx in range(-nimg,nimg+1):
        for ny in range(-nimg,nimg+1):
            R=D+np.array([nx*Wx,ny*Wy]); r=np.hypot(R[:,0],R[:,1])
            with np.errstate(divide="ignore",invalid="ignore"): t=np.where(r>1e-12,erfc(alpha*r)/r,0.0)
            out+=t
    zero=np.hypot(D[:,0],D[:,1])<1e-12; out[zero]-=2*alpha/np.sqrt(np.pi)
    kmax=10*alpha; mx=int(np.ceil(kmax*Wx/(2*np.pi))); my=int(np.ceil(kmax*Wy/(2*np.pi)))
    for p in range(-mx,mx+1):
        for q in range(-my,my+1):
            if p==0 and q==0: continue
            k=np.array([2*np.pi*p/Wx,2*np.pi*q/Wy]); kn=np.hypot(*k)
            if kn>kmax: continue
            out+=(2*np.pi/(A*kn))*erfc(kn/(2*alpha))*np.cos(D@k)
    out-=2*np.sqrt(np.pi)/(A*alpha)
    return out
def phi_matrix(X,I,J,W,Lx,Ly,**kw):
    """Full N x N Ewald matrix using translation invariance (2 reference sites: even row and odd row)."""
    N=len(X); P=np.zeros((N,N))
    for ref in (0,Lx):    # site (0,0) and (0,1)
        D=X-X[ref]; D[:,0]-=np.round(D[:,0]/W[0])*W[0]; D[:,1]-=np.round(D[:,1]/W[1])*W[1]
        ph=ewald_phi(D,W,**kw)
        rows=np.where(J%2==(ref//Lx))[0]
        for i in rows:
            di=I[i]; dj=J[i]-J[ref]
            jj=((J+dj)%Ly)*Lx+((I+di)%Lx)      # site that the reference's partner maps to
            P[i,jj]=ph
    return P
@numba.njit(cache=True)
def kawasaki_nn(m,NB,T,nsteps,seed):
    np.random.seed(seed); N=len(m)
    for t in range(len(T)):
        beta=1.0/T[t]
        for s in range(nsteps):
            i=np.random.randint(N); j=np.random.randint(N)
            if m[i]==m[j]: continue
            hi=0.0; hj=0.0; adj=0
            for k in range(6):
                hi+=m[NB[i,k]]; hj+=m[NB[j,k]]
                if NB[i,k]==j: adj=1
            dmi=m[j]-m[i]; dmj=-dmi
            dE=-(dmi*(hi-adj*m[j])+dmj*(hj-adj*m[i]))
            if dE<=0 or np.random.random()<np.exp(-beta*dE):
                tmp=m[i]; m[i]=m[j]; m[j]=tmp
    return m
@numba.njit(cache=True)
def kawasaki_full(m,P,T,nsteps,seed):
    np.random.seed(seed); N=len(m); h=P@m
    for t in range(len(T)):
        beta=1.0/T[t]
        for s in range(nsteps):
            i=np.random.randint(N); j=np.random.randint(N)
            if m[i]==m[j]: continue
            dmi=m[j]-m[i]; dmj=-dmi
            # U = -1/2 sum_{ij} m_i m_j P_ij (diag = self images); exact change
            dE=-(dmi*h[i]+dmj*h[j]+0.5*dmi*dmi*P[i,i]+0.5*dmj*dmj*P[j,j]+dmi*dmj*P[i,j])
            if dE<=0 or np.random.random()<np.exp(-beta*dE):
                tmp=m[i]; m[i]=m[j]; m[j]=tmp
                for k in range(N): h[k]+=P[k,i]*dmi+P[k,j]*dmj
    return m
def U_nn(m,E): return -np.sum(m[E[:,0]]*m[E[:,1]])
def U_full(m,P): return -0.5*m@P@m
