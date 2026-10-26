# rev7 regeneration of rev6 Fig 16 (CoE density modulation m(r)). Same shape law and coefficients as the ray trace actually used
# in the rev7 solve (famP.py 'poly', C=[-0.041,-0.035,0.1353,0], run by joint_solve/phase2/p2_bg_build.py at D=34.96):
# m(r)=exp sum_k c_k (min(r/D,R_b)-1)^(k+1), R_b=1.5(1+R_V/D), R_V=14.30048 Gpc (observable radius). rev6 drew it at D=42.
import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
CC=[-0.041,-0.035,0.1353,0.0];RV=14.30048287260589
def law(D):
    Rb=1.5*(1+RV/D);return (lambda r:np.exp(sum(c*(np.minimum(np.abs(r/D),Rb)-1)**(k+1) for k,c in enumerate(CC)))),Rb
fig,ax=plt.subplots(figsize=(13,4.8));r=np.linspace(0.01,80,1600);out={}
for D,st,c in ((34.96,'-','b'),):
    m,Rb=law(D);s=(m(D+1e-3)-m(D-1e-3))/2e-3;re=min(D+RV,Rb*D)
    ri=r[r<=re];ax.plot(ri,m(ri),st,color=c,lw=2.2 if D<40 else 1.2,label=f'D={D:g} Gpc'+f': m(D)={m(D):.4f}, slope at observer {s:.2e}/Gpc, R_b D={Rb*D:.1f} Gpc')
    ax.plot(D,m(D),'o',color=c,ms=8);out[str(D)]=dict(m_obs=float(m(D)),slope=float(s),m_lo=float(m(D-RV)),m_hi=float(m(re)),RbD=float(Rb*D))
ax.axvspan(34.96-RV,34.96+RV,color='gold',alpha=.3,label=f'rev7 observable range D±R_V = {34.96-RV:.1f}–{34.96+RV:.1f} Gpc')
ax.axvspan(34.96+RV,80,color='0.88',alpha=.6,label='beyond observable edge: unconstrained (no continuation drawn)')
ax.set_xlim(0,80);ax.set_xlabel('distance from CoE [Gpc]');ax.set_ylabel('m(r) (before normalising)');ax.legend(fontsize=8,loc='lower left');ax.grid(alpha=.3)
ax.set_title('CoE density modulation m(r) on the rev7 background (cut at 80 Gpc)',fontsize=10)
fig.tight_layout();fig.savefig(__file__.replace('.py','.png'),dpi=140,facecolor='w');json.dump(out,open(__file__.replace('.py','.json'),'w'),indent=1);print(out)
