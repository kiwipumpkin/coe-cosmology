"""Void size function and cluster mass function from the same linear sigma(R) used elsewhere in clump_scan
(density_vs_scale/lcdm_sigma_tophat.txt: EH98 no-wiggle, s8=0.81, Om=0.31, h=0.70; R in Mpc, no h).
  void size function : Jennings, Li & Hu 2013 (arXiv:1304.6087) volume-conserving 'Vdn' model on top of the
                       Sheth & van de Weygaert 2004 (astro-ph/0311260) two-barrier excursion set:
                       delta_v = -2.7 (linear), delta_c = 1.686, r = 1.7 r_L (rho_v = 0.2 rho_m)
  cluster mass function: Tinker et al. 2008 (arXiv:0803.2706), Delta = 200 x mean, z = 0:
                       A=0.186, a=1.47, b=2.57, c=1.19
All abundances per Mpc^3 (no h)."""
import os,numpy as np
CS="/workspace/clump_scan"
_T=np.loadtxt(os.path.join(CS,"density_vs_scale","lcdm_sigma_tophat.txt"))
LR=np.log(_T[:,0]); LS=np.log(_T[:,1])
RHO_M=0.31*2.775e11*0.70**2            # Msun / Mpc^3
MSUN=1.989e30
def sigma(R):
    x=np.log(np.asarray(R,float)); sl=(LS[-1]-LS[-2])/(LR[-1]-LR[-2]); s0=(LS[1]-LS[0])/(LR[1]-LR[0])
    return np.exp(np.where(x>LR[-1],LS[-1]+sl*(x-LR[-1]),np.where(x<LR[0],LS[0]+s0*(x-LR[0]),np.interp(x,LR,LS))))
def dlnsinv_dlnR(R,e=1e-3): return -(np.log(sigma(R*(1+e)))-np.log(sigma(R*(1-e))))/(np.log(1+e)-np.log(1-e))
DV=-2.7; DC=1.686; RFAC=1.7
def f_svdw(s,nterm=200):
    """SvdW two-barrier multiplicity f_lnsigma (their eq. via Jennings+13 eq. 6)."""
    D=abs(DV)/(DC+abs(DV)); x=D*np.atleast_1d(s)/abs(DV); j=np.arange(1,nterm+1)[:,None]
    return np.sum(2*j*np.pi*x**2*np.sin(j*np.pi*D)*np.exp(-(j*np.pi*x)**2/2),0)
def void_dn_dlnr(r):
    """Vdn: dn/dlnr = f(sigma(r_L))/V(r) * dln(1/sigma)/dln r_L,  r_L = r/1.7"""
    r=np.atleast_1d(np.asarray(r,float)); rL=r/RFAC; V=4/3*np.pi*r**3
    return f_svdw(sigma(rL))/V*dlnsinv_dlnR(rL)
def void_n_above(rmin,rmax=1e5):
    lr=np.linspace(np.log(rmin),np.log(rmax),4000); y=void_dn_dlnr(np.exp(lr)); return float(np.trapezoid(y,lr))
def void_volfrac_above(rmin,rmax=1e5):
    lr=np.linspace(np.log(rmin),np.log(rmax),4000); r=np.exp(lr); y=void_dn_dlnr(r)*4/3*np.pi*r**3; return float(np.trapezoid(y,lr))
TA,Ta,Tb,Tc=0.186,1.47,2.57,1.19
def R_of_M(M): return (3*np.asarray(M,float)/(4*np.pi*RHO_M))**(1/3)
def cluster_dn_dlnM(M):
    """Tinker08 200m: dn/dlnM = f(sigma) rho_m/M dln(1/sigma)/dlnM"""
    M=np.atleast_1d(np.asarray(M,float)); s=sigma(R_of_M(M)); f=TA*((s/Tb)**(-Ta)+1)*np.exp(-Tc/s**2)
    return f*RHO_M/M*dlnsinv_dlnR(R_of_M(M))/3.0
def cluster_n_above(Mmin,Mmax=1e18):
    lm=np.linspace(np.log(Mmin),np.log(Mmax),3000); return float(np.trapezoid(cluster_dn_dlnM(np.exp(lm)),lm))
def sample_void_radii(rng,n,rmin,rmax=None):
    rmax=rmax or rmin*30; lr=np.linspace(np.log(rmin),np.log(rmax),3000); p=void_dn_dlnr(np.exp(lr)); c=np.cumsum(p); c/=c[-1]
    return np.exp(np.interp(rng.uniform(size=n),c,lr))
def cluster_masses_ranked(n,Mmin,V):
    """abundance matching: the k-th most massive of n clusters in volume V has n(>M)V = k - 1/2"""
    lm=np.linspace(np.log(Mmin),np.log(1e17),3000); M=np.exp(lm)
    N=np.array([cluster_n_above(m) for m in M[::30]])*V; lmS=lm[::30]
    k=np.arange(1,n+1)-0.5
    return np.exp(np.interp(-np.log(k),-np.log(np.maximum(N,1e-30)),lmS))
if __name__=="__main__":
    h=0.7
    for rm in (10,15,20,30,40,60,80):
        print(f"voids r>{rm} Mpc ({rm*h:.0f} Mpc/h): n={void_n_above(rm):.3e}/Mpc^3, vol frac {void_volfrac_above(rm):.3f}")
    for Mm in (1e14,3e14,1e15,3e15):
        print(f"clusters M200m>{Mm:.0e} Msun: n={cluster_n_above(Mm):.3e}/Mpc^3")
    V=4/3*np.pi*14300.0**3; print("visible sphere volume Mpc^3",f"{V:.3e}")
