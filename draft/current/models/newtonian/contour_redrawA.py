import numpy as np
src=open("contour_h70.py").read(); exec(src.split("# Map A")[0])
d=np.load("contour_A.npz"); RD,ms,yA,Z,bad=d["RD"],d["ms"],d["yA"],d["Z"],d["bad"]
keep=(yA>=0)&(yA<=3.2); ms,yA,Z,bad=ms[keep],yA[keep],Z[keep],bad[keep]
X,Y=np.meshgrid(RD,yA); dots=[(0.246,np.interp(m,ms,yA)) for m in (1,4,10)]
draw(X,Y,Z,bad,"R / D  (galaxy distance / CoE distance)","v_esc / v  at observer","RD_vesc",dots,["1×","4×","10×"])
ok=~bad
def thr(i,cond):
    out=[]
    for j in range(len(RD)):
        s=np.where(ok[:,j]&cond[:,j])[0]; out.append(yA[s].min() if len(s) else np.nan)
    return np.array(out)
for nm,i,cond in (("sd>10",1,Z[...,1]>10),("mean<65",0,Z[...,0]<65)):
    t=thr(i,cond); print(nm,[(round(r,2),round(x,2)) for r,x in zip(RD[::5],t[::5])])
print("max valid R/D per row:",[(round(y,2),RD[ok[i]].max() if ok[i].any() else None) for i,y in enumerate(yA)][::4])
B=np.load("contour_B.npz"); mB,kB,ZB,bB=B["mB"],B["kB"],B["Z"],B["bad"]
for a in range(0,len(kB),3): print("k=%.2f valid m<=%s mean range %.0f-%.0f sd max %.0f"%(kB[a],mB[~bB[a]].max() if (~bB[a]).any() else None,ZB[a][~bB[a],0].min() if (~bB[a]).any() else 0,ZB[a][~bB[a],0].max() if (~bB[a]).any() else 0,ZB[a][~bB[a],1].max() if (~bB[a]).any() else 0))
a=np.argmin(abs(kB-1)); print("k=1: m where sd>10:",mB[(ZB[a,:,1]>10)&~bB[a]].min(),"mean<65:",mB[(ZB[a,:,0]<65)&~bB[a]].min() if ((ZB[a,:,0]<65)&~bB[a]).any() else None)
