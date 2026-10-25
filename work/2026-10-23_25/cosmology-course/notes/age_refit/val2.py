import numpy as np,ltb2 as T,forward as F,sys,json
sys.path.insert(0,'/workspace/coe_grz'); import grav_redshift as GRZ
out={}
def obs(P):
    L=T.LTB(P,'fwd'); o=F.integrate(P,True); io=o['io']; return L,o,np.log(F.R0[io]),o['t']*F.C
# 1 FLRW limit: uniform density (exp law, L=1e8 Gpc >> all labels), negligible point mass (1e40 kg), p=1
P=F.params_from('exp',40,0.6,1.0,13,8); L,o,xo,tn=obs(P); r=[]
for th in np.radians([0,45,90,135,180]):
    s=T.ray(L,tn,xo,th,0.05*tn,rtol=1e-11); z=np.exp(s.y[2,-1])-1
    a=L.R(s.t[-1],np.array([xo]))[0]/L.R(tn,np.array([xo]))[0]; zf=1/a-1
    # analytic FLRW comoving distance check: chi = int dt/a  -> compare with emitter label distance (x_e - x_o along radial for th=0)
    r.append(dict(theta=float(np.degrees(th)),z=z,z_FLRW=zf,relerr=(z-zf)/zf,status=s.status))
    print('FLRW',r[-1],flush=True)
out['FLRW']=r
# 2 Schwarzschild: point mass dominant, negligible dust (f=1e-12); radial rays vs Cosmology's exact one_plus_z
P=F.params_from('exp',53,0.8,1.0,-12,8); L,o,xo,tn=obs(P); r=[]
rs=lambda x: 2*L.m(x)
for th,sk in ((0.0,-1),(np.pi,1)):
    s=T.ray(L,tn,xo,th,(0.3 if th==0 else 0.9)*tn,rtol=1e-11); xe=s.y[0,-1]; te=s.t[-1]; z=np.exp(s.y[2,-1])-1
    Re=L.R(te,np.array([xe]))[0]; Ue=L.Rdot(te,np.array([xe]))[0]; Ro=L.R(tn,np.array([xo]))[0]; Uo=L.Rdot(tn,np.array([xo]))[0]
    RS=2*float(L.m(np.log(1e-3))); ee=np.sqrt(Ue**2+1-RS/Re); eo=np.sqrt(Uo**2+1-RS/Ro)
    zs=GRZ.one_plus_z(RS,(Re,ee,0.0,1),(Ro,eo,0.0,1),0.0,sk,sk)-1 if s.ncross==0 else np.nan-1
    r.append(dict(theta=float(np.degrees(th)),z=z,z_schw=zs,relerr=(z-zs)/zs,rs_over_Ro=RS/Ro,Ue=Ue,Uo=Uo)); print('SCHW',r[-1],flush=True)
out['Schwarzschild_radial']=r
json.dump(out,open('gr2/validation.json','w'),indent=1,default=float)
