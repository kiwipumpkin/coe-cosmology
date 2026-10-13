"""Simple4 model level sets of signed LOS relative velocity (km/s) vs angle (0 = away from CoE).
v(r)^2 = v0^2 + 2GM(1/r - 1/D), sign of v0; invalid (grey) where v^2<0."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
G=6.674e-11; GPC=3.0857e25; M0=7.5e52; D0=58.2; R0=14.3
vesc=lambda M,D:np.sqrt(2*G*M/(D*GPC))
ang=np.linspace(0,360,721); t=np.radians(ang); ct,st=np.cos(t),np.sin(t)
def vl(D,R,M,v0):
    X=D+R*ct; Y=R*st; r=np.hypot(X,Y); q=v0**2+2*G*M*(1/(r*GPC)-1/(D*GPC))
    v=np.where(q>=0,np.sign(v0)*np.sqrt(np.abs(q)),np.nan); return (v*(X*ct+Y*st)/r-v0*ct)/1e3
VB=0.4*vesc(M0,D0); print("baseline |v0| = %.0f km/s"%(VB/1e3))
def zeros(Z,ys,yv):
    i=np.argmin(abs(ys-yv)); z=Z[i]; ok=~np.isnan(z); j=np.where(np.diff(np.sign(z[ok])))[0]; return np.round(ang[ok][j],1)
def plot(ys,Z,ylab,title,fn,logy=False,mark=None):
    f,ax=plt.subplots(figsize=(9,5.2)); ax.set_facecolor("0.85")
    ax.pcolormesh(ang,ys,np.where(np.isnan(Z),np.nan,0),cmap="Greys",vmin=0,vmax=10,shading="auto")
    hi,lo=np.nanmax(Z),np.nanmin(Z)
    if hi>0: cp=ax.contour(ang,ys,Z,levels=np.linspace(0,hi,8)[1:-1],colors="#c0392b",linewidths=1); ax.clabel(cp,fmt="%.0f",fontsize=8)
    if lo<0: cn=ax.contour(ang,ys,Z,levels=np.linspace(lo,0,8)[1:-1],colors="#2e5fa8",linewidths=1); ax.clabel(cn,fmt="%.0f",fontsize=8)
    c0=ax.contour(ang,ys,Z,levels=[0],colors="k",linewidths=2.6); ax.clabel(c0,fmt="0",fontsize=10)
    if mark is not None: ax.axhline(mark,color="k",lw=.8,ls=":")
    if logy: ax.set_yscale("log")
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE",fontsize=11)
    ax.set_ylabel(ylab,fontsize=11); ax.set_title(title+"\nsigned LOS relative velocity (km/s): red receding, blue approaching, bold = 0, grey = v² < 0",fontsize=10.5)
    f.savefig(fn,dpi=130,bbox_inches="tight"); plt.close(f)
# (1) signed x
xs=np.linspace(-3,3,241); xs=xs[abs(xs)>0.05]
Z=np.array([vl(D0,R0,M0,x*vesc(M0,D0)) for x in xs])
plot(xs,Z,"velocity field: our shell speed / escape speed\n(all shells same energy; − collapsing, + expanding)","D = 58.2, R = 14.3 Gpc, M = 7.5×10⁵² kg (v_esc = 74,700 km/s)","ls_vel_angle_v0.png",mark=0)
for x in (-3,-1,-0.6,-0.4,0.4,0.6,1,3): print("x=",x,"zeros",zeros(Z,xs,x))
for sg,nm in ((-1,"col"),(1,"exp")):
    v0=sg*VB; w="collapsing" if sg<0 else "expanding"
    Rs=np.linspace(1,30,117); Z=np.array([vl(D0,R,M0,v0) for R in Rs])
    plot(Rs,Z,"galaxy distance R (Gpc)",f"{w}, v0 = {sg*VB/1e3:,.0f} km/s (x = 0.4; sets whole velocity field, all shells same energy), D = 58.2, M = 7.5×10⁵² kg",f"ls_vel_angle_R_{nm}.png",mark=14.3)
    print(nm,"R zeros:",{R:zeros(Z,Rs,R).tolist() for R in (1,7,14.3,25)})
    Ds=np.linspace(20,80,121); Z=np.array([vl(D,R0,M0,v0) for D in Ds])
    plot(Ds,Z,"our distance to CoE, D (Gpc)",f"{w}, v0 = {sg*VB/1e3:,.0f} km/s fixed (sets whole velocity field), R = 14.3, M = 7.5×10⁵² kg",f"ls_vel_angle_D_{nm}.png",mark=58.2)
    print(nm,"D zeros:",{D:zeros(Z,Ds,D).tolist() for D in (20,40,58.2,80)})
    Ms=np.logspace(-1,1,121); Z=np.array([vl(D0,R0,M0*k,v0) for k in Ms])
    plot(Ms,Z,"M / 7.5×10⁵² kg  (v0 fixed in km/s)",f"{w}, v0 = {sg*VB/1e3:,.0f} km/s fixed (sets whole velocity field), D = 58.2, R = 14.3",f"ls_vel_angle_M_{nm}.png",logy=True,mark=1)
    print(nm,"M zeros:",{k:zeros(Z,Ms,k).tolist() for k in (0.1,0.5,1,2,10)})
