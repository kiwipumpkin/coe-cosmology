"""Linear LCDM sigma(R) at z=0: Eisenstein & Hu (1998) no-wiggle transfer, Omega_m=0.31, Omega_b=0.049, h=0.70, n_s=0.965, sigma8=0.81.
Gaussian (sigma_G) and top-hat (sigma_TH) smoothing. R in Gpc (physical, h=0.7 -> converted to Mpc/h)."""
import numpy as np
from scipy.integrate import simpson
Om,Ob,h,ns,s8=0.31,0.049,0.70,0.965,0.81; Tcmb=2.7255
def T_nw(k):          # k in h/Mpc
    omh2=Om*h*h; fb=Ob/Om; th=Tcmb/2.7
    s=44.5*np.log(9.83/omh2)/np.sqrt(1+10*(Ob*h*h)**0.75)
    aG=1-0.328*np.log(431*omh2)*fb+0.38*np.log(22.3*omh2)*fb**2
    q=k*th**2/(Om*h); kh=k*h
    G=Om*h*(aG+(1-aG)/(1+(0.43*kh*s)**4))
    q=k*th**2/G
    L=np.log(2*np.e+1.8*q); C=14.2+731/(1+62.5*q); return L/(L+C*q*q)
k=np.logspace(-6,2,20000); P0=k**ns*T_nw(k)**2
def W_th(x): return 3*(np.sin(x)-x*np.cos(x))/x**3
def _sig(Rh,W): return np.sqrt(simpson(P0*W(k*Rh)**2*k**2/(2*np.pi**2),x=k))
A=(s8/_sig(8.0,W_th))**2
def sigma_TH(R_gpc): return np.sqrt(A)*_sig(R_gpc*1000*h,W_th)
def sigma_G(R_gpc): return np.sqrt(A)*_sig(R_gpc*1000*h,lambda x:np.exp(-x*x/2))
if __name__=="__main__":
    print("check sigma8 TH(8/h Mpc)=",sigma_TH(8/1000/h))
    for R in [0.3,0.5,1,2,4]: print(f"R={R} Gpc  sigma_Gauss={sigma_G(R):.2e}  sigma_tophat={sigma_TH(R):.2e}")
