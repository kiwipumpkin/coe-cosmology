import numpy as np,warnings;warnings.filterwarnings('ignore')
from scipy.integrate import quad;from scipy.optimize import brentq
RV=14.30048287260589;OM=0.3;DEL=0.031
def shape(kind,x,p):
    if kind=='exp':return np.exp(-x)
    if kind=='eta':return np.exp(-np.abs(x)**p)
    if kind=='soft':return 1/(1+np.abs(x)**p)
def Lock(kind,D,p):   # L (Gpc) with [rho(D-d)-rho(D+d)]/rho(D)=0.031, d=14.3 Gpc
    f=lambda lL:(shape(kind,(D-RV)/np.exp(lL),p)-shape(kind,(D+RV)/np.exp(lL),p))/shape(kind,D/np.exp(lL),p)-DEL
    g=np.linspace(-3,12,400);v=np.array([f(x) for x in g]);ok=np.isfinite(v);s=np.where(ok[:-1]&ok[1:]&(np.sign(v[:-1])!=np.sign(v[1:])))[0]
    if len(s)==0:return np.nan
    i=s[-1];return np.exp(brentq(f,g[i],g[i+1]))   # largest L root (gentlest gradient)
def age(R,M,E):return quad(lambda r:1/np.sqrt(2*M/r+2*E),0,R,limit=200)[0]
def ratio(kind,D,p):   # units D=1, uD=1 (cancels); common bang; returns H_par/H_perp, M(<D)/Mcrit, L
    L=Lock(kind,D,p)
    if not np.isfinite(L):return np.nan,np.nan,np.nan
    lD=L/D;s1=shape(kind,1/lD,p)
    M=lambda R:1.5*OM*quad(lambda r:shape(kind,r/lD,p)/s1*r*r,0,R,limit=200)[0]
    M1=M(1.);t0=age(1.,M1,0.5-M1)
    def u(R):
        m=M(R);E=brentq(lambda E:age(R,m,E)-t0,-m/R*0.999999,60);return np.sqrt(2*m/R+2*E)
    h=1e-3;return (u(1+h)-u(1-h))/(2*h),M1*2,L
