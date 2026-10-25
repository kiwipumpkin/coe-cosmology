import numpy as np
def add(u,v):  # u (+) v : velocity v measured in frame moving at u
    u=np.asarray(u,float);v=np.asarray(v,float);uu=u@u
    if uu==0: return v
    g=1/np.sqrt(1-uu); up=(u@v)/uu*u; vp=v-(v@u)/uu*u
    return (u+ up*0 + (v@u)/uu*u - (v@u)/uu*u + up)/(1+u@v) + vp/(g*(1+u@v)) if False else (u+up)/(1+u@v)+vp/(g*(1+u@v))
def rest_grad(Hpar,Hperp,bf,D,eps=1e-6):
    bx=np.array([bf,0.]); g=1/np.sqrt(1-bf**2)
    def w(x):
        r=np.array([x[0]+D,x[1]]); rh=r/np.linalg.norm(r)
        vs=np.array([Hpar*x[0],Hperp*x[1]]); return add(vs,bf*rh)
    def up(x):  # velocity in origin-member rest frame and transformed position
        u=w(x); ux,uy=u; d=1-bf*ux
        return np.array([(ux-bf)/d, uy/(g*d)]), np.array([g*x[0],x[1]])
    G=np.zeros((2,2))
    for j in range(2):
        e=np.zeros(2);e[j]=eps; u1,x1=up(e); u0,x0=up(-e)
        G[:,j]=(u1-u0)/(x1-x0)[j]
    return G
if __name__=="__main__":
    D=1.0
    for bf in [0,0.3,0.6,0.9]:
        G=rest_grad(1.0,0.9,bf,D); g=1/np.sqrt(1-bf**2)
        print(bf,G.round(6),"gamma",round(g,4),"bf/D/g",round(bf/D/g,4))
