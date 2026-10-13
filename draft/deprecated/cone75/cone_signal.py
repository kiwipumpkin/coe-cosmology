import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
src=open("cone_aaall.py").read().split("sty=")[0].replace('Ks={"raw":None,','Ks={')
exec(src)
AALL=1.114  # Migkas+2020 Eq.6
SIG={"60":(0.051/AALL,0.069/AALL),"75":(0.050/AALL,0.061/AALL)}   # Sect 5.3.3 / 5.3.2 extreme-cone sigmas
for mf in (1,4):
  for c in ("60","75"):
    f,ax=plt.subplots(figsize=(9,5)); lo,hi=SIG[c]
    ax.axhspan(1-hi,1+hi,color="0.88",zorder=0,label="Migkas 2020 per-cone σ(A)/A_all, %s° cone: ±%.3f–%.3f"%(c,lo,hi))
    ax.axhspan(1-lo,1+lo,color="0.78",zorder=0)
    ys=[]
    for sg,col,nm in ((1,"#c0392b","expanding"),(-1,"#2e5fa8","collapsing")):
        for k,ls in (("%s° top-hat"%c,"-"),("%s° cos taper"%c,"--")):
            h=full(Ks[k]@A(sg,mf)[:i180]); ys.append(h); ax.plot(ang,h,color=col,ls=ls,lw=2,label="%s, %s"%(nm,k))
    for k in ("top-hat","cos taper"):
        e=full(Ks["%s° %s"%(c,k)]@A(1,mf)[:i180]); q=full(Ks["%s° %s"%(c,k)]@A(-1,mf)[:i180]); d=e-q
        print("M=%g %s° %-9s diff 0°=%+.3f 90°=%+.3f 180°=%+.3f | in σ(lo..hi): 0° %.1f..%.1f  90° %.1f..%.1f  180° %.1f..%.1f"%(mf,c,k,d[0],d[360],d[720],abs(d[0])/hi,abs(d[0])/lo,abs(d[360])/hi,abs(d[360])/lo,abs(d[720])/hi,abs(d[720])/lo))
    a=np.concatenate(ys); pad=.15*(a.max()-a.min()); ax.set_ylim(min(a.min()-pad,1-hi-.005),max(a.max()+pad,1+hi+.005))
    ax.axhline(1,color="0.4",lw=.8)
    ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel("angle θ (deg):  0° = away from CoE, 180° = toward CoE")
    ax.set_ylabel("cone-averaged A / A_all"); ax.set_title("Cone-averaged A/A_all, %s° cones, point mass %g× (D = 58.2, R = 14.3 Gpc)\ngrey: Migkas+2020 A&A 636 A15 per-cone σ (Sect. 5.3.%s), dark = smaller, light = larger"%(c,mf,"3" if c=="60" else "2"),fontsize=10)
    ax.legend(frameon=False,fontsize=8.5,ncol=2,loc="upper center",bbox_to_anchor=(0.5,-0.13))
    f.savefig("cone_signal_M%g_%s.png"%(mf,c),dpi=120,bbox_inches="tight"); plt.close(f)
