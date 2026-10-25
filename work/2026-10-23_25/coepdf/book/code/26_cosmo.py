import numpy as np
h=0.7;Om=0.31;Ob=0.049;ns=0.965;S8=0.8;dc=1.686;dv=-2.7
rho_m=Om*2.775e11*h**2   # Msun/Mpc^3
def T_EH(k):  # Eisenstein-Hu 1998 no-wiggle, k in 1/Mpc
    om=Om*h*h;ob=Ob*h*h;fb=Ob/Om;s=44.5*np.log(9.83/om)/np.sqrt(1+10*ob**0.75)
    a=1-0.328*np.log(431*om)*fb+0.38*np.log(22.3*om)*fb**2;th=2.7255/2.7
    G=Om*h*(a+(1-a)/(1+(0.43*k*s)**4));q=k*th**2/G/h*h  # q=k/(h Mpc^-1)*Theta^2/Gamma ; k already 1/Mpc -> k/h
    q=k/h*th**2/G
    L=np.log(2*np.e+1.8*q);C=14.2+731/(1+62.5*q);return L/(L+C*q*q)
def W(x):return 3*(np.sin(x)-x*np.cos(x))/x**3
_k=np.geomspace(1e-5,50,6000)
def _sig(R,A=1.0):return np.sqrt(np.trapezoid(A*_k**(ns+2)*T_EH(_k)**2*W(_k*R)**2/(2*np.pi**2),_k))
A=(S8/_sig(8/h))**2
def P(k):return A*k**ns*T_EH(k)**2  # Mpc^3
def sigR(R):return _sig(R,A)
def MofR(R):return 4/3*np.pi*R**3*rho_m
def RofM(M):return (3*M/(4*np.pi*rho_m))**(1/3)
def ps_dndlnM(M):
    s=np.array([sigR(RofM(m)) for m in M]);ls=np.gradient(np.log(s),np.log(M))
    return np.sqrt(2/np.pi)*rho_m/M*(dc/s)*np.abs(ls)*np.exp(-dc**2/(2*s*s)),s
def svdw_dndlnR(Rnl):  # Sheth & van de Weygaert 2004, linear-theory void size function (dv=-2.7), R_nl=1.7 R_L
    RL=Rnl/1.7;s=np.array([sigR(r) for r in RL]);D=abs(dv)/(dc+abs(dv));x=D*s/abs(dv)
    j=np.arange(1,200)[:,None];f=2*np.sum(j*np.pi*x**2*np.sin(j*np.pi*D)*np.exp(-(j*np.pi*x)**2/2),0)
    dls=np.gradient(np.log(1/s),np.log(RL));return f/(4/3*np.pi*RL**3)*dls   # per lnR_L = per lnR_nl
if __name__=='__main__':
    print('A',A,'sig8',sigR(8/h),'sig(100/h)',sigR(100/h),'sig(140/h)',sigR(140/h))
    M=np.geomspace(1e14,2e15,30);n,s=ps_dndlnM(M);print('n(>1e14) Mpc^-3',np.trapezoid(n,np.log(M)),'frac mass',np.trapezoid(n*M,np.log(M))/rho_m)
    R=np.geomspace(10,50,30);nv=svdw_dndlnR(R);print('nvoid Mpc^-3',np.trapezoid(nv,np.log(R)),'vol frac',np.trapezoid(nv*4/3*np.pi*R**3,np.log(R)))
    kk=np.geomspace(1e-4,1,400);print('Pk peak k',kk[np.argmax(P(kk))],'-> scale 2pi/k Mpc',2*np.pi/kk[np.argmax(P(kk))])
