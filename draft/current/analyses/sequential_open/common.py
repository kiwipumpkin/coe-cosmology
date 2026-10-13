import os,sys,pickle,json,numpy as np
SEQ="/workspace/clump_scan/sequential"; OUT=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,SEQ); sys.path.insert(0,OUT)
_c=os.getcwd(); os.chdir(SEQ)
import seq_core
from seq_core import C_G,A_IN_G,R,FLOOR,rms,amps,RMIN
os.chdir(_c)
import open_region as O
O.setup(C_G,A_IN_G)
V_OLD_B=268960.0; V_NEW_B=635768.2      # 14.3-42.9 minus core  vs  ball(42.9) U core
V_OLD_A=88275.2;  V_NEW_A=455900.3      # 14.3-30   minus core  vs  ball(30)   U core
SUF={"both":"","core":"_core","inside":"_inside","tidal":"_tidal","tidalcore":"_tidalcore"}[O.MODE]
V_REG={"both":635768.2,"core":635768.2-12143.0,"inside":12143.0+268787.0,"tidal":None,"tidalcore":None}[O.MODE]
def patch(mod):
    mod.region_ok=O.region_ok; mod.rand_in_shell=O.rand_in_region
def save_atomic(obj,fn):
    tmp=fn+".tmp"
    with open(tmp,"wb") as f: pickle.dump(obj,f)
    os.replace(tmp,fn)
def where(X,r=None):
    X=np.atleast_2d(X); s=np.linalg.norm(X,axis=1); dc=np.linalg.norm(X-C_G,axis=1)
    return dict(n=int(len(X)),inside_us=int((s<RMIN).sum()),in_core=int((dc<A_IN_G).sum()),beyond_42_9=int((s>42.9).sum()))

if O.MODE in ("tidal","tidalcore"):    # ball(30) minus the CoE core (tidal) / whole ball(30) (tidalcore)
    _r=np.random.default_rng(0); _P=_r.uniform(-30,30,(3_000_000,3)); _s=np.linalg.norm(_P,axis=1); _d=np.linalg.norm(_P-C_G,axis=1)
    V_REG=float(60**3*np.mean((_s<=30)&((_d>A_IN_G) if O.MODE=="tidal" else True)))
    V_CORE_IN=float(60**3*np.mean((_s<=30)&(_d<=A_IN_G)))

def frame30(ax,lw=0.8,fs=9):
    """Default side-view frame (Chris, 10 Oct): only the modeled region, our sphere + shell to 30 Gpc;
    the CoE direction is a VECTOR from us (origin) out to the 30 Gpc edge, labelled 'to CoE' near its tip."""
    from matplotlib.patches import Circle as _C
    ax.add_patch(_C((0,0),30.0,fill=False,ec="0.4",lw=lw,ls=":",zorder=4))
    u=C_G[:2]/np.linalg.norm(C_G[:2]); w=np.array([-u[1],u[0]])
    ax.annotate("",xy=(30.0*u[0],30.0*u[1]),xytext=(0,0),arrowprops=dict(arrowstyle="-|>,head_length=0.8,head_width=0.4",color="#a93226",lw=2.0,shrinkA=0,shrinkB=0),zorder=8)
    p=22.0*u+2.6*w
    ax.text(p[0],p[1],"to CoE",color="#a93226",fontsize=fs,fontweight="bold",ha="center",va="center",zorder=8,bbox=dict(fc="white",ec="none",alpha=0.8,pad=0.4),gid="coe_label")
    ax.set_xlim(-32,32); ax.set_ylim(-32,32); ax.set_aspect("equal")
