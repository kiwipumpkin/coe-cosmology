import sys; sys.path.insert(0,'.'); from grcore import *
mu=0.1; v0=0.99; to=106.0
oat,_=geodesic(mu,np.zeros(2),np.array([v0,0]),0.0,to); Xo,Uo=oat(to); ro,po=polar(Xo); wo=static_vel(mu,ro,po,Uo)
for d in (0,60,92.5,130,180):
    t=np.radians(d); e=np.array([np.cos(t),np.sin(t)]); at,_=geodesic(mu,e,relsum(v0,0,e),0.0,to)
    r,_=observe(mu,Xo,Uo,[at],to,0.0)
    # recover solution by re-solving and extracting
    g=lambda s: np.linalg.norm(Xo-at(s)[0])-(to-s); s0=brentq(g,1e-9,to-1e-6)
    Xe=at(s0)[0]; dd=Xe-Xo; dr=dd@np.array([np.cos(po),np.sin(po)]); dp=dd@np.array([-np.sin(po),np.cos(po)])
    def res(z):
        so,L=ray(mu,Xo,z[0],to,0.0)
        if so.t[-1]>z[1]: return [5,5]
        rr,pp,_=so.sol(z[1]); return cart(rr,pp)-at(z[1])[0]
    b,s=fsolve(res,[np.arctan2(-dp,-dr/np.sqrt(1-2*mu/ro)),s0],xtol=1e-12)
    so,L=ray(mu,Xo,b,to,0.0); re,pe,pr=so.sol(s); fe=1-2*mu/re; fo=1-2*mu/ro
    Xe,Ue=at(s); we=static_vel(mu,re,pe,Ue)
    # photon direction in local static frames (unit vectors, cartesian axes)
    def ndir(r,p,pr,L):
        f=1-2*mu/r; kr=pr/np.sqrt(f)*1; kp=L/r*np.sqrt(f)/1  # static-frame components (scaled)
        n=np.array([kr,kp])/np.hypot(kr,kp); return np.array([n[0]*np.cos(p)-n[1]*np.sin(p), n[0]*np.sin(p)+n[1]*np.cos(p)])
    ne=ndir(re,pe,pr,L); no=ndir(ro,po,np.cos(b),L)
    De=1/np.sqrt(1-we@we)*(1-we@ne); Do=1/np.sqrt(1-wo@wo)*(1-wo@no)
    bend=np.degrees(np.arccos(np.clip(ne@no,-1,1)))
    print(f"{d:6.1f}deg emit t={s:.3f} at x=({Xe[0]-0:.2f},{Xe[1]:.2f}) dist-from-mass r={re:.3f} | emitter speed {np.hypot(*we):.5f} obs speed {np.hypot(*wo):.5f} | angle emitter-vel to photon {np.degrees(np.arccos(we@ne/np.hypot(*we))):.3f}deg | light bent {bend:.4f}deg | grav sqrt(fe/fo)={np.sqrt(fe/fo):.4f} Doppler obs/emit={Do/De:.4f} total={r[0]:.4f} check={np.sqrt(fe/fo)*Do/De:.4f}",flush=True)
