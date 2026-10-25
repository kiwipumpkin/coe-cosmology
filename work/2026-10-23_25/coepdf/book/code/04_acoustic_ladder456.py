# Ladder R=4,5,6 Gpc with per-R exact/Limber correction (13 nodes l 51-400, log-ratio interpolated). Field only.
import os,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.special import spherical_jn
import acoustic_model as M
from acoustic_model import AFwd,Gc,rho,T0,LMAX
import acoustic_fit as A
from acoustic_ladder import lim,pk,MULT,p0,PlP,pl1
CK='acoustic_ladder_456.json';S=json.load(open(CK)) if os.path.exists(CK) else {};ll=np.arange(LMAX+1)
NODES=(51,60,75,90,110,130,150,175,200,240,280,330,400)
def exact(F,p,L):
    ch,KM,dch=F.ch,F.KM,F.dch;W=np.exp(-0.5*((ch-p['chiw'])/p['sw'])**2);Wd=-W*(ch-p['chiw'])/p['sw']**2;m=W>1e-4;ch,KM,W,Wd=ch[m],KM[m],W[m],Wd[m]
    k=np.linspace(L/ch.max()*0.7,min((L+60)/ch.min()*2.5,4000),2500);x=k[:,None]*ch[None]
    Aa=spherical_jn(L,x)@(KM*W)*dch;Bp=spherical_jn(L+1,x)@(KM*Wd)*dch;Bm=spherical_jn(L-1,x)@(KM*Wd)*dch
    acc=k**2*Aa**2/3+(Bp**2+Bm**2)/6+k*Aa*(Bp-Bm)/3;Pp=(4*np.pi*Gc*rho)**2*F.Pd(p,k)/k**4
    return 2/np.pi*np.trapezoid(k**2*Pp*acc,k)*L*(L+1)/(2*np.pi)*T0**2
for R in (4.0,5.0,6.0):
    key=f'R{R:g}'
    if key in S and 'D1' in S[key]:continue
    F=AFwd(R=R,lmax=M.LX,kmax=200,nk=500);p=dict(p0,chiw=R-0.056,sw=0.05);t=F.terms(p)
    Dl=lim(F,p,np.arange(M.LX+1,LMAX+1));D1=np.r_[0,0,t['scalar']+t['mlens'],Dl]
    rr=S.get(key,{}).get('ratio',[])
    for L in NODES[len(rr):]:
        rr.append(float(exact(F,p,L)/lim(F,p,np.array([L]))[0]));S[key]=dict(ratio=rr);json.dump(S,open(CK,'w'))   # checkpoint per node
    c=np.ones(LMAX+1);sel=(ll>M.LX)&(ll<=400);c[sel]=np.exp(np.interp(ll[sel],NODES,np.log(rr)));D1=D1*c
    h=float(D1[150:350].max());lev=MULT*h;i=[j for j in range(8) if lev[j]<=pl1<=lev[j+1]]
    S[key]=dict(ratio=rr,D1=D1.tolist(),peaks=pk(D1),first_peak_height_x1=h,first_peak_l=int(150+np.argmax(D1[150:350])),bracket=[(float(MULT[j]),float(MULT[j+1])) for j in i],mult_to_match=pl1/h,jump_l50=float(D1[51]/D1[50]))
    json.dump(S,open(CK,'w'));print(key,S[key]['peaks'],'h %.4g'%h,S[key]['bracket'],'mult %.3g'%(pl1/h),'D51/D50 %.3f'%S[key]['jump_l50'],np.round(rr,3).tolist(),flush=True)
fig,axs=plt.subplots(1,3,figsize=(21,7));cs=plt.cm.viridis(np.linspace(0,1,9));l=ll[2:]
for ax,R in zip(axs,(4,5,6)):
    s=S[f'R{R:g}'];D1=np.array(s['D1'])
    ax.plot(l,A.Lt[2:],'k-',lw=1.8,label='ΛCDM');ax.plot(A.Pf[:,0],A.Pf[:,1],'.',color='grey',ms=1.5,alpha=.6,label='Planck 2018')
    for j,m in enumerate(MULT):
        Dm=m*D1[2:];Dr=Dm*np.random.default_rng(42+j).chisquare(2*l+1)/(2*l+1);ax.plot(l,Dm,'-',color=cs[j],lw=.6,alpha=.3);ax.plot(l,Dr,'-',color=cs[j],lw=.6,label=f'×{m:.3g}')
    ax.set_xscale('log');ax.set_yscale('log');ax.set_ylim(1,1e6);ax.set_xlim(2,2500);ax.set_xlabel('ℓ');ax.set_ylabel('D_ℓ (μK²)')
    ax.set_title(f'R={R} Gpc, shell {R-0.056:.3f}±0.05; peaks {s["peaks"]}\n1st-peak multiplier {s["mult_to_match"]:.2f} (bracket ×{s["bracket"][0][0]:.3g}–×{s["bracket"][0][1]:.3g})',fontsize=9);ax.legend(loc='lower center',ncol=6,fontsize=6,framealpha=.85)
fig.suptitle(f'acoustic free_W ladder, exact/Limber-corrected (one realization per level, faint = mean). Planck peaks {PlP}',fontsize=10);fig.tight_layout();fig.savefig('acoustic_ladder_456.png',dpi=100)
