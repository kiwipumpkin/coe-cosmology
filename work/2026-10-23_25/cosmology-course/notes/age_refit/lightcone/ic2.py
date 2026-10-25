"""Intercept construction v2 (Chris 2026-10-16): global shell motion U>c (expansion of the shell, a(t)=1+U t, d1=1, c=1) carries
observer and emitters; PECULIAR velocities relative to the shell (< c) enter D_motion; D_grav as before.
Photon moves at c relative to the local shell -> straight line in comoving coords with comoving path tau=(1/U) ln(1+U t).
Peculiar motion (constant v_pec, radial from O) also moves comoving distance v_pec*tau, so the intercept equation is the old one in tau:
  (v_o^2-1) tau^2 + 2 dx v_o tau + |d|^2 = 0   ->  t_f = (exp(U tau)-1)/U,   global redshift factor 1+z_U = a(t_f) = 1+U t_f.
S = D_motion * D_grav / (1+U t_f).   Peculiar speeds: v_o^2=beta1^2+2mu(1/s-1/s)=beta1^2; v_e^2 = beta1^2 re0^2 + 2 mu (1/re - 1/re0)  (old law, now PECULIAR).
Rejected: v_pec^2 outside (0,1), 2mu/r>=1, U<=1."""
import numpy as np
def S_ic2(beta1,q,mu,s,U,ell_rel):
    th=np.radians(ell_rel); xe0=1+q*np.cos(th); ye0=q*np.sin(th); re0=np.hypot(xe0,ye0)
    if U<=1 or np.any(re0<1e-6): return None
    ro,re=s,s*re0; vo2=beta1**2; ve2=beta1**2*re0**2+2*mu*(1/re-1/re0)
    if not(0<vo2<1) or np.any(ve2<=0) or np.any(ve2>=1) or 2*mu/ro>=1 or np.any(2*mu/re>=1): return None
    vo,ve=np.sqrt(vo2),np.sqrt(ve2); dx,dy=1-xe0,-ye0; a=vo2-1; b=2*dx*vo; c=dx**2+dy**2
    tau=(-b-np.sqrt(b*b-4*a*c))/(2*a); nx,ny=(dx+vo*tau)/tau,dy/tau; bex,bey=ve*xe0/re0,ve*ye0/re0
    Dm=np.sqrt(1-ve2)/np.sqrt(1-vo2)*(1-vo*nx)/(1-(bex*nx+bey*ny)); Dg=np.sqrt((1-2*mu/re)/(1-2*mu/ro))
    tf=np.expm1(U*tau)/U; return dict(S=Dm*Dg/(1+U*tf),Dm=Dm,Dg=Dg,zU=U*tf,tau=tau,tf=tf,xe0=xe0,ye0=ye0,vo=vo,ve=ve)
