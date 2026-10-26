import numpy as np, json
from scipy import optimize, integrate, stats
p=print
aeq=1/3400
F=lambda a:(2/3)*(a-2*aeq)*np.sqrt(a+aeq)+(4/3)*aeq**1.5
num=integrate.quad(lambda x:x/np.sqrt(aeq+x),0,0.5)[0]; p("eq2 closed vs num",F(0.5),num)
Hs=F(1)/13.8e9
p("t_eq kyr",F(aeq)/Hs/1e3,"t_rec kyr",F(1/1090)/Hs/1e3)
p("implied H0 (km/s/Mpc) EdS-like: H*=",Hs*977.8e9, "H0=H*sqrt(1+aeq)",Hs*977.8e9*np.sqrt(1+aeq))
p("EdS age at H0=70 Gyr",2/3/(70/977.8))
c=299792.458
for H0 in [70,47]:
  chi=2*c/H0*(1-1/np.sqrt(1091))/1000; p("matter-only comoving dist to LSS Gpc, H0",H0,chi)
p("H0 needed for chi*=4 Gpc", 2*c*(1-1/np.sqrt(1091))/4000)
p("V_R",4/3*np.pi*64,"V_D+R",4/3*np.pi*62.2**3,"ratio%",64/62.2**3*100,"Vout",4/3*np.pi*72.83**3,"r%",64/72.83**3*100)
for l in [2,10,30,220,1000]: p("CV",l,np.sqrt(2/(2*l+1)))
p("Dl ex",10*11/2/np.pi*1e-10*(2.7255e6)**2)
d=14.3005
for dl in [0.031,2*0.01554]:
  L=d/np.arcsinh(dl/2); p("L exp",dl,L,0.3*np.exp(58.2/L),-1/L)
p("L d=3",3/np.arcsinh(0.0155))
laws={'power':lambda x:(1+x)**-1.2,'cored':lambda x:1/(1+x*x),'gauss':lambda x:np.exp(-x*x/2),
'soft':lambda x:(1+x**0.3)**-1,'stretch':lambda x:np.exp(-x**0.25),'etaedge':lambda x:np.exp(-x**0.05)}
paperL={'power':1049.3,'cored':321.9,'gauss':231.5,'soft':5109.6,'stretch':15246,'etaedge':0.849}
D=58.2
for k,g in laws.items():
  f=lambda L:(g((D-d)/L)-g((D+d)/L))/g(D/L)-0.031
  L=paperL[k]; dlt=f(L)+0.031
  h=1e-4; grad=(np.log(g((D+h)/L))-np.log(g((D-h)/L)))/(2*h)
  p(k,"paperL",L,"delta(L)",dlt,"rho0",0.3/g(D/L),"dlnrho/dr",grad)
G=4.30e-12; rho=0.31*2.775e11*0.49*1e9; p("rho",rho,"4piGrho",4*np.pi*G*rho,"Phi/c2",4*np.pi*G*rho*1e-4/c**2)
p("Md void",0.85*4.215e10*4/3*np.pi*30**3)
p("P ex",10**-8.282*100**-0.183*np.exp(-(100/335)**2),"w/o damping",10**-8.282*100**-0.183)
p("dl",2*np.pi*3.944/0.0819,2*np.pi*3.5715/0.0795,300*0.0795/2/np.pi,[x*np.pi/220*1000 for x in (4,8,14)])
p("r2",2*np.pi*3.94/593*1000,"kB",990/3.94,"Rc/w",3.5715/0.3355,(4-3.5715)/0.3355)
x=56.3/42-1; cc=[-0.041,-0.035,0.1353,0]; p("m(56.3)",np.exp(sum(ci*x**(i+1) for i,ci in enumerate(cc))),"slope",cc[0]/42, "Rb",1.5*(1+1/2.937))
for dof in [2,6,12]: p("chi2 levels 2p",stats.chi2.ppf([.6827,.9545,.9973],2))
# Migkas
m=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True)
y=m['A_over_Aall']; s=(m['A_err_low']+m['A_err_high'])/2
for nm,sig in [('mean',s),('low',m['A_err_low']),('high',m['A_err_high']),('side',np.where(y>1,m['A_err_low'],m['A_err_high']))]:
  p("null chi2",nm,((y-1)**2/sig**2).sum())
p("N bins",len(y),"spacing",np.diff(m['l_deg'])[:3],"p(chi2>45.14|71)",stats.chi2.sf(45.14,71),"p(7.67|69)",stats.chi2.cdf(7.67,69),"p(10.8|68) cdf",stats.chi2.cdf(10.8,68))
p("chi2/dof",7.67/(71-3),10.79/(71-4),11.6/(71-2))
# effective independent bins: 90deg cones step 5 -> ~360/90=4
Neff=4; p("binomial 11/12",stats.binom.sf(10,12,.5))
p("geometric residual after 95.9% cancel of 1-5%", [x*0.041 for x in (0.01,0.05)])
p("Planck chi2 dof", 3787/(2499-12), 2558/(2499-6),"dchi2",3787-2558)
