import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
exec(open("rel_vel_collapse_combo_c.py").read().split("sets=")[0])
vesc=lambda M,D: np.sqrt(2*GM0*M/(D*GPC))/c
def calc(M,D,R,v):
    ve=vesc(M,D)
    if not (0<=v<ve): return np.nan,np.nan,1
    if R>=0.95*D: return np.nan,np.nan,2
    y=v/ve; rst=np.inf if y>=1 else D/(1-y*y); yy=curve(M,D,R,rst)
    return 100*(w*(yy>0)).sum()/w.sum(),(w*yy).sum()/w.sum(),0
D0,M0,R0,V0=58.2,1.0,14.3,0.1
ve0=vesc(1,D0)
# name, xlabel, ylabel, xr, yr, map(x,y)->(M,D,R,v), held
P=[("RD_vvesc","R / D","|v| / v_esc(D)",(0.02,0.9),(0.01,0.99),lambda x,y:(M0,D0,x*D0,y*vesc(M0,D0)),"M=1×, D=58.2 Gpc"),
   ("RD_stallD","R / D","stall radius / D",(0.02,0.9),(1.01,20),lambda x,y:(M0,D0,x*D0,np.sqrt(1-1/y)*ve0),"M=1×, D=58.2 Gpc"),
   ("RD_energy","R / D","v² D / M  (c² Gpc per 7.5×10⁵² kg)",(0.02,0.9),(0.01,3.5),lambda x,y:(M0,D0,x*D0,np.sqrt(y*M0/D0)),"M=1×, D=58.2 Gpc (v varies)"),
   ("DM_Rv","D / M  (Gpc per 7.5×10⁵² kg)","R / |v|  (Gpc per c)",(10,600),(10,500),lambda x,y:(M0,x*M0,y*V0,V0),"M=1×, |v|=0.1c"),
   ("MD_RD","M / D  (7.5×10⁵² kg per Gpc)","R / D",(0.003,0.3),(0.02,0.9),lambda x,y:(x*D0,D0,y*D0,V0),"D=58.2 Gpc, |v|=0.1c"),
   ("Mv2_RD","M / v²  (7.5×10⁵² kg per c²)","R / D",(10,2000),(0.02,0.9),lambda x,y:(M0,D0,y*D0,np.sqrt(M0/x)),"M=1×, D=58.2 Gpc (v varies)"),
   ("Dv_Rv","D / |v|  (Gpc per c)","R / |v|  (Gpc per c)",(100,6000),(10,800),lambda x,y:(M0,x*V0,y*V0,V0),"M=1×, |v|=0.1c"),
   ("DM_RM","D / M  (Gpc per 7.5×10⁵² kg)","R / M  (Gpc per 7.5×10⁵² kg)",(10,600),(1,80),lambda x,y:(M0,x*M0,y*M0,V0),"M=1×, |v|=0.1c"),
   ("DM_Rtimesv","D / M  (Gpc per 7.5×10⁵² kg)","R · |v|  (Gpc · c)",(10,600),(0.1,8),lambda x,y:(M0,x*M0,R0,y/R0),"M=1×, R=14.3 Gpc (v varies)"),
   ("vvesc_stallD_check","|v| (c)","M  (× 7.5×10⁵² kg)",(0.005,0.5),(0.1,10),lambda x,y:(y,D0,R0,x),"D=58.2, R=14.3 Gpc")]
N=42; pdf=PdfPages("heat_combo_gallery.pdf"); summ=[]
def panel(ax,X,Y,Z,code,qty,lg,xl,yl,title):
    for k,g in ((1,"0.6"),(2,"0.85")): ax.pcolormesh(X,Y,np.where(code==k,1,np.nan),cmap=matplotlib.colors.ListedColormap([g]),shading="auto")
    if qty=="recede": pc=ax.pcolormesh(X,Y,Z,vmin=0,vmax=100,cmap="RdBu_r",shading="auto"); lev=50; cl="% sky receding"
    else:
        a=np.nanmax(abs(Z)); pc=ax.pcolormesh(X,Y,Z,norm=matplotlib.colors.TwoSlopeNorm(0,-max(a,1e-9),max(a,1e-9)) if np.nanmin(Z)<0<np.nanmax(Z) else matplotlib.colors.TwoSlopeNorm(0,-a,a),cmap="RdBu_r",shading="auto"); lev=0; cl="sky-average v_LOS / c"
    nseg=0
    if np.nanmin(Z)<lev<np.nanmax(Z):
        cs=ax.contour(X,Y,Z,[lev],colors="k",linewidths=1.4); nseg=len(cs.allsegs[0])
    if lg: ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(xl,fontsize=8); ax.set_ylabel(yl,fontsize=8); ax.tick_params(labelsize=7)
    cb=plt.colorbar(pc,ax=ax); cb.set_label(cl,fontsize=7); cb.ax.tick_params(labelsize=6); ax.set_title(title,fontsize=8)
    return nseg
for name,xl,yl,xr,yr,mp,held in P:
    G=plt.figure(figsize=(11,8.5)); axs=G.subplots(2,2); res={}
    for j,lg in enumerate((False,True)):
        sp=np.geomspace if lg else np.linspace; X=sp(*xr,N); Y=sp(*yr,N)
        F=np.full((N,N),np.nan); S=F.copy(); C=np.zeros((N,N))
        for a,yv in enumerate(Y):
            for b,xv in enumerate(X): F[a,b],S[a,b],C[a,b]=calc(*mp(xv,yv))
        for i,(qty,Z) in enumerate((("recede",F),("skymean",S))):
            sc="log" if lg else "lin"; t="%s, %s axes  [held: %s]"%(("% receding" if qty=="recede" else "sky mean v/c"),sc,held)
            ns=panel(axs[i,j],X,Y,Z,C,qty,lg,xl,yl,t)
            f1,a1=plt.subplots(figsize=(6.5,5)); panel(a1,X,Y,Z,C,qty,lg,xl,yl,"Collapse (stall-radius free fall): "+t); f1.savefig("heat_combo/heat_combo_%s_%s_%s.png"%(name,qty,sc),dpi=110,bbox_inches="tight"); plt.close(f1)
            res[(qty,sc)]=(np.nanmin(Z),np.nanmax(Z),ns)
    G.suptitle("Pair: %s   (grey dark: v ≥ v_esc; light: R ≥ 0.95 D)"%name,fontsize=11); G.tight_layout(); pdf.savefig(G); plt.close(G)
    print(name,{k:("%.3g..%.3g seg%d"%v) for k,v in res.items()})
pdf.close()
