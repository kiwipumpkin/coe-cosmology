import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.optimize import brentq
exec(open("wr.py").read().split("out={}")[0].replace('plt.rcParams.update','_=dict'))
plt.rcParams.update({"font.size":16,"axes.labelsize":16,"xtick.labelsize":14,"ytick.labelsize":14})
SEED=0.01; SPAN=45; mid=lambda t,a: a
pages=[]
for m in (1,10):
    M=m*M1; r0=brentq(lambda x:run(x,M,800)[1][-1]-D,0.05*D,D); t,r,v=run(r0,M); g=G*M/r**2
    for fr in (0.016,0.03):
        lam=3*fr*D; ell=2*np.pi*D/lam            # fixed angular mode
        kk=2*np.pi/lam*D/r; Ncum=np.concatenate([[0],np.cumsum(0.5*(np.sqrt(kk*g)[1:]+np.sqrt(kk*g)[:-1])*np.diff(t))])
        amp=np.minimum(SEED*np.exp(Ncum),1.0)      # in units of Delta, capped at Delta
        idx=[0,len(t)//2,len(t)-1]; cols=["#9ecae1","#4292c6","#08306b"]
        W,HW,HG=8,5.0,1.9; fig=plt.figure(figsize=(W,HW+2*HG+2.4))
        H=HW+2*HG+2.4; aw=fig.add_axes([0.04,(2*HG+1.2)/H,0.92,HW/H]); aw.axis("off"); Rm=D*1.04/GPC; aw.set_xlim(-0.05,Rm+0.05); aw.set_ylim(-0.25,Rm*np.sin(np.radians(SPAN))+0.05)
        th=np.radians(np.linspace(0,SPAN,3000))
        aw.plot([0,Rm],[0,0],color="0.6",lw=1); aw.plot([0,Rm*np.cos(th[-1])],[0,Rm*np.sin(th[-1])],color="0.6",lw=1)
        for i,c in zip(idx,cols):
            R_=r[i]; Dl=fr*R_; dr=amp[i]*Dl*np.cos(ell*th)
            ro=(R_+0.5*Dl+dr)/GPC; ri=(R_-0.5*Dl+dr)/GPC
            aw.fill(np.r_[ro*np.cos(th),(ri*np.cos(th))[::-1]],np.r_[ro*np.sin(th),(ri*np.sin(th))[::-1]],color=c,lw=0)
            aw.text(R_/GPC*np.cos(np.radians(47)),R_/GPC*np.sin(np.radians(47)),"%.1f Gyr"%(t[i]/GYR),fontsize=14,ha="right",color=c)
        aw.plot(0,0,"k*",ms=16); aw.text(0.08,-0.15,"CoE",fontsize=14)
        fig.text(0.5,1-0.75/H,"%d× mass,  Δ/R = %.3f,  λ = 3Δ = %.0f Mpc now (ℓ ≈ %.0f)\nseed ripple %.2gΔ, grown by RT e-folds, capped at Δ"%(m,fr,lam/MPC,ell,SEED),fontsize=14,va="center",ha="center")
        aw.text(D/GPC*0.55,-0.18,"radius (Gpc): shell at 0 → %.2f, mid → %.2f, now → %.2f"%(r[0]/GPC,r[idx[1]]/GPC,D/GPC),fontsize=12,ha="center")
        fig.text(0.5,1-0.25/H,"Our shell wrinkling (RT growth, Vishniac spacing)",ha="center",fontsize=17,weight="bold")
        a1=fig.add_axes([0.16,(HG+1.0)/H,0.78,(HG-0.35)/H]); deg=np.degrees(th)
        for i,c in zip(idx,cols): a1.plot(deg,amp[i]*np.cos(ell*th),color=c,lw=2.5)
        a1.set_xlim(0,SPAN); a1.set_ylim(-1.1,1.1); a1.set_ylabel("ripple\nδr / Δ"); a1.axhline(0,color="0.6",lw=.7); a1.set_xlabel("angle along shell (deg)")
        a2=fig.add_axes([0.16,0.55/H,0.78,(HG-0.35)/H])
        a2.plot(t/GYR,amp,color="#08306b",lw=3); a2.set_yscale("log"); a2.set_ylim(SEED*0.8,1.5); a2.set_ylabel("amplitude\n/ Δ")
        for i,c in zip(idx,cols): a2.plot(t[i]/GYR,amp[i],"o",color=c,ms=10)
        a2.set_xlabel("time (Gyr)"); a2.text(0.02,0.85,"e-folds by now: %.1f"%Ncum[-1],transform=a2.transAxes,fontsize=14)
        fn="wrinkle_shell_M%d_d%03d.png"%(m,round(fr*1000)); fig.savefig(fn,dpi=130); plt.close(fig); pages.append(fn)
        print(fn,"Nnow %.2f amp %s"%(Ncum[-1],[round(amp[i],3) for i in idx]))
from PIL import Image
ims=[Image.open(p).convert("RGB") for p in pages]; ims[0].save("wrinkle_shell.pdf",save_all=True,append_images=ims[1:],resolution=130)
