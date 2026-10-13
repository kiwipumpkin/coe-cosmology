import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
exec(open("levelD_v2_nomass.py").read().split("Lv=")[0])   # run(sg,D,mf) with sheet method, k, R=14.3
D=58.2; Ms=[0,0.5,1,2,4]; cols=plt.cm.viridis(np.linspace(0.1,0.9,4))
mean=lambda y:(w*y).sum()/w.sum(); rec=lambda y:100*(w*(y>0)).sum()/w.sum()
res={(s,m):run(s,D,m) for s in (1,-1) for m in Ms}
nm={1:"expanding (V = +k r)",-1:"collapsing (V = −k r)"}
for idx,tag,yl in ((0,"vlos","v_LOS / c  (+ = receding)"),(1,"aaall","A / A_all = |±kR n + dt Δg| / kR")):
    lo=min(res[x][idx].min() for x in res); hi=max(res[x][idx].max() for x in res); pad=.05*(hi-lo)
    for s,sn in ((1,"expand"),(-1,"collapse")):
        f,ax=plt.subplots(figsize=(9,5.3)); ax.axhline(0 if idx==0 else 1,color="0.6",lw=.8)
        for m in Ms:
            y=res[s,m][idx]; lab="M = %g× (strength %.3g) | sky mean %+.4f%s"%(m,0.247*m,mean(y)," c, rec %.0f%%"%rec(y) if idx==0 else "")
            if m==0: ax.plot(ang,y,"k--",lw=1.6,label="no mass | "+lab.split("| ")[1])
            else: ax.plot(ang,y,color=cols[Ms.index(m)-1],lw=2,label=lab)
            print(tag,sn,m,"min %.4f max %.4f mean %.4f"%(y.min(),y.max(),mean(y)))
        ax.set_ylim(lo-pad,hi+pad); ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE"); ax.set_ylabel(yl)
        ax.set_title("%s — %s, mass level sets\nsheet method (aaall.py), D = 58.2, R = 14.3 Gpc; M in units of 7.5×10⁵² kg"%(nm[s],"LOS velocity" if idx==0 else "A/A_all"),fontsize=10)
        ax.legend(frameon=False,fontsize=8.5,loc="upper center",bbox_to_anchor=(.5,-.13),ncol=1)
        f.savefig("mass_%s_%s_c.png"%(tag,sn),dpi=125,bbox_inches="tight"); plt.close(f)
