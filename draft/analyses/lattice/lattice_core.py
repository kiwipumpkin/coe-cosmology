"""Lattice cluster/void structure in the surround region (same region, density law, tidal code and residual as
../shell_ensemble/surround and ../shell_ensemble/with_clumps). Nothing upstream is modified; code is imported by exec."""
import os, json, numpy as np
HERE=os.path.dirname(os.path.abspath(__file__)); SE=os.path.join(os.path.dirname(HERE),"shell_ensemble")
_cwd=os.getcwd()
g={"__file__":os.path.join(SE,"with_clumps","ensemble_with_clumps.py"),"__name__":"imp"}
os.chdir(os.path.join(SE,"with_clumps")); exec(open(g["__file__"]).read().split('\nif __name__')[0],g); os.chdir(_cwd)
field3,amps,R,rms,RV,FLOOR,KG,radius_gpc,C,baseline,P4,kk=[g[k] for k in
   ("field3","amps","R","rms","RV","FLOOR","KG","radius_gpc","C","baseline","P4","kk")]
sg=g["g"]; A_IN,B_SH,S_MAX,L,D=[sg[k] for k in ("A_IN","B_SH","S_MAX","L","D")]
GPC=3.0857e25; RHO=2.9e-27          # same constants as ensemble.py
BASE0=rms(R)/FLOOR                   # 237.3
R4,_sub4=baseline(4); BASE4=rms(R4)/FLOOR   # 88.8 (four clumps, compensated)
AX=C/np.linalg.norm(C)               # unit vector observer -> CoE (lattice z-axis = CoE axis)

# ---- unit cells: fractional basis positions and signs (+1 cluster, -1 void); 'a' = cubic/hex cell edge
CELLS={
 "SC_rocksalt": dict(vec=np.eye(3),basis=np.array([[0,0,0]]),sign=None,desc="simple cubic, neighbours alternate +/- (rock-salt / NaCl sign pattern); a = +/- neighbour spacing"),
 "BCC_CsCl":   dict(vec=np.eye(3),basis=np.array([[0,0,0],[.5,.5,.5]]),sign=np.array([1,-1]),desc="body-centred cubic: clusters at corners, voids at centres (CsCl); sign of amplitude free = reverse; a = cube edge"),
 "FCC_L10":    dict(vec=np.eye(3),basis=np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]]),sign=np.array([1,1,-1,-1]),desc="face-centred cubic, alternating (001) layers +/- (CuAu L1_0); a = cube edge"),
 "HCP_AB":     dict(vec=np.array([[1,0,0],[.5,np.sqrt(3)/2,0],[0,0,np.sqrt(8/3)]]),basis=np.array([[0,0,0],[1/3,1/3,.5]]),sign=np.array([1,-1]),desc="hexagonal close-packed, A layers clusters / B layers voids; a = in-plane spacing, c = 1.633a"),
}
def _rot(e):
    a,b,c=e; ca,sa,cb,sb,cc,sc=np.cos(a),np.sin(a),np.cos(b),np.sin(b),np.cos(c),np.sin(c)
    Rz1=np.array([[ca,-sa,0],[sa,ca,0],[0,0,1]]); Ry=np.array([[cb,0,sb],[0,1,0],[-sb,0,cb]]); Rz2=np.array([[cc,-sc,0],[sc,cc,0],[0,0,1]])
    return Rz1@Ry@Rz2
# frame whose z-axis is the CoE axis (Euler angles are relative to it)
_z=AX; _x=np.cross([0,0,1.0],_z); _x/=np.linalg.norm(_x); _y=np.cross(_z,_x); FRAME=np.c_[_x,_y,_z]

def in_region(P):
    s=np.linalg.norm(P,axis=1); return (s>1.0)&(s<=S_MAX)&(np.linalg.norm(P-C,axis=1)>=A_IN)

def _frac_points(kind,a,Rm,phase):
    """Fractional lattice points (n+basis+phase) whose positions can lie within S_MAX of the observer."""
    cell=CELLS[kind]; V=cell["vec"]*a
    o=np.linalg.solve(V.T,Rm.T@(-C))-np.asarray(phase)          # observer in lattice (cell-vector) coordinates
    span=int(np.ceil(S_MAX/(a*min(np.linalg.norm(cell["vec"],axis=1))*0.5)))+2   # generous (non-orthogonal cells)
    r=[np.arange(int(np.floor(o[i]))-span,int(np.floor(o[i]))+span+1) for i in range(3)]
    n=np.stack(np.meshgrid(*r,indexing="ij"),-1).reshape(-1,3)
    pts=[];sg_=[]
    for j,b in enumerate(cell["basis"]):
        pts.append(n+b+phase)
        sg_.append(np.where(n.sum(1)%2==0,1,-1) if cell["sign"] is None else np.full(len(n),cell["sign"][j]))
    return np.concatenate(pts),np.concatenate(sg_),V

def lattice_sites(kind,a,euler,phase,clip=True):
    """a in R_vis units; euler (3, rad, ZYZ) relative to the CoE-axis frame; phase (3) in cell fractions. Lattice origin at the CoE (+phase)."""
    if kind=="shell_rings": return shell_rings(a,*phase,clip=clip)
    Rm=FRAME@_rot(euler)
    f,s,V=_frac_points(kind,a,Rm,np.asarray(phase,float))
    P=C+(f@V)@Rm.T
    if not clip: return P,s
    ok=in_region(P); return P[ok],s[ok]

def shell_rings(a,u_r,u_th,u_ph,clip=True):
    """Shell-adapted lattice about the CoE axis: concentric spheres around the CoE (radial spacing a), on each sphere rings at
    polar angle steps ~a/r about the CoE->observer axis, ring members evenly spaced (arc ~a), signs alternate along ring, ring-to-ring
    and sphere-to-sphere. Free: a, radial offset u_r, polar offset u_th, azimuthal rotation u_ph (all fractions of a step)."""
    pts=[];sg_=[]; j0=int(np.floor((A_IN-u_r*a)/a))
    for j in range(j0, j0+int((D+S_MAX-A_IN)/a)+3):
        r=(j+u_r)*a
        if r<A_IN or r>D+S_MAX: continue
        nth=max(1,int(round(np.pi*r/a))); dth=np.pi/nth
        # observer sits at polar angle pi (CoE->observer axis points from C toward origin); only rings with polar angle near pi can be in region
        for k in range(nth):
            t=(k+0.5+(u_th-0.5))*dth
            if t<=0 or t>=np.pi: continue
            n=max(1,int(round(2*np.pi*r*np.sin(t)/a)))
            if n%2: n+=1                                   # even so +/- alternation closes around the ring
            ph=(np.arange(n)+u_ph)*2*np.pi/n
            loc=np.c_[r*np.sin(t)*np.cos(ph),r*np.sin(t)*np.sin(ph),-r*np.cos(t)*np.ones(n)]   # local z = observer->CoE; -cos so t~pi is near observer
            sgn=np.where(np.arange(n)%2==0,1,-1)*(1 if (j+k)%2==0 else -1)
            pts.append(loc); sg_.append(sgn)
    P=np.concatenate(pts); s=np.concatenate(sg_)
    P=C+P@FRAME.T*(-1)*np.array([1,1,1])   # flip so that local -z (polar pi) points from C toward observer
    if not clip:
        k=np.linalg.norm(P,axis=1)<=S_MAX+2*a; return P[k],s[k]
    ok=in_region(P); return P[ok],s[ok]

def masses(P,s,a,kind,amp):
    """Mass per object = amp * rho(r_CoE) * V_site (fraction of the cell's own mass moved), clusters +, voids -;
    voids then rescaled so total void deficit = total cluster excess (net zero, same convention as the random runs).
    Returns M (kg), void_scale, V_site (m^3)."""
    nb=1 if kind in ("SC_rocksalt","shell_rings") else len(CELLS[kind]["basis"])
    Vsite=(a*RV*GPC)**3*(abs(np.linalg.det(CELLS[kind]["vec"])) if kind!="shell_rings" else 1.0)/nb
    rho=RHO*np.exp(-(np.linalg.norm(P-C,axis=1)-D)/L)          # exp-lock, normalised to mean density at D
    base=rho*Vsite
    Mp=base*(s>0); Mn=base*(s<0)
    vs=Mp.sum()/Mn.sum() if Mn.sum()>0 else 1.0
    return amp*(Mp-vs*Mn), vs, Vsite

def pattern(P,s,a,kind):
    """Unit-amplitude tidal pattern (amp=1) and bookkeeping."""
    M,vs,Vs=masses(P,s,a,kind,1.0)
    if len(P)==0: return np.zeros_like(R),vs,Vs,M
    Dd=np.linalg.norm(P,axis=1)
    w=np.zeros_like(R)
    for i in range(0,len(P),4000):
        sl=slice(i,i+4000); w+=(M[sl]/KG)@field3(np.arctan2(P[sl,1],P[sl,0]),np.arcsin(np.clip(P[sl,2]/Dd[sl],-1,1)),Dd[sl])
    return w,vs,Vs,M

def amp_max(vs): return 1.0/max(1.0,vs)       # cluster excess <= cell mass and void deficit <= cell mass (delta >= -1)

def best_amp(target,w,vs):
    ww=float(np.dot(w,w))
    if ww==0: return 0.0,rms(target)
    A=float(np.dot(target,w)/ww); Am=amp_max(vs); A=float(np.clip(A,-Am,Am))
    return A,rms(target-A*w)
