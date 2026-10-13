from ltbx import *; import sys
D=float(sys.argv[1]); rs=float(sys.argv[2]); M=Model(D,rs,rc=0.1,NQ=6000,NT=3000)
w=1/np.sqrt(1+M.E2)
# exact radial inward-looking (180deg) null ray, backward: dr/dtau = -c/X
ts=np.concatenate([np.linspace(M.t[-1],M.t[0],200000)])
r=M.ro; rr=[r]; tt=[ts[0]]
for k in range(len(ts)-1):
    dt=ts[k]-ts[k+1]; _,_,Rp,_=M.f(ts[k],np.array([r])); X=Rp[0]*np.interp(r,M.r,w)
    r=r-C*dt/X; 
    if r<=M.r[1]: break
    rr.append(r); tt.append(ts[k+1])
rr=np.array(rr); tt=np.array(tt)
for Rt in [float(x) for x in sys.argv[3].split(',')]:
    # find emitter label where present proper length = Rt
    X0=np.interp(rr,M.r,M.Rp[-1])*np.interp(rr,M.r,w); L=np.concatenate([[0],np.cumsum(-np.diff(rr)*0.5*(X0[1:]+X0[:-1]))])
    if L[-1]<Rt: print("D=%g rs=%g R=%g: beyond reach (max %.2f)"%(D,rs,Rt,L[-1])); continue
    j=np.searchsorted(L,Rt); seg=slice(0,j+1); rs_=rr[seg]; dr=-np.diff(rs_); mid=lambda a:.5*(a[1:]+a[:-1])
    Xd_g=np.array([M.f(t,np.array([q]))[3][0] for t,q in zip(tt[seg],rs_)])*np.interp(rs_,M.r,w)
    Xd_0=np.interp(rs_,M.r,M.Rdp[-1])*np.interp(rs_,M.r,w); X0s=X0[seg]
    Rp0=np.interp(rs_,M.r,M.Rp[-1]); Rdp0=np.interp(rs_,M.r,M.Rdp[-1])
    num_g=(mid(Xd_g)*dr).sum(); num_0=(mid(Xd_0)*dr).sum(); den=(mid(X0s)*dr).sum()
    Hc=num_g/den*KM; HA=num_0/den*KM; Hcl=(mid(Rdp0)*dr).sum()/(mid(Rp0)*dr).sum()*KM
    print("D=%g rs=%g R=%g: H_cosmo=%.1f | same but neighbours' separation speed taken NOW=%.1f | classical (areal, Euclid)=%.1f | t_emit=%.2f Gyr, emitter CoE areal now=%.2f Gpc; mean ratio Rdot'(t_pass)/Rdot'(now)=%.3f"%(D,rs,Rt,Hc,HA,Hcl,tt[j],np.interp(rs_[-1],M.r,M.R[-1]),num_g/num_0),flush=True)
