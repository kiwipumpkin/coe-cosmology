import numpy as np, healpy as hp
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':15})
ns=32;L=2
def quad(alm_dict):
 a=np.zeros(hp.Alm.getsize(L),complex)
 for m,v in alm_dict.items():a[hp.Alm.getidx(L,2,m)]=v
 return hp.alm2map(a,ns)
Q1=quad({0:1.0,1:0.3+0.2j,2:0.4-0.1j})
Q2=quad({0:0.9,1:0.25+0.25j,2:0.3-0.2j})   # similar orientation
al=0.8;S=Q1-al*Q2
fig=plt.figure(figsize=(15,8));v=np.abs(Q1).max()
for i,(m,t) in enumerate([(Q1,'Q₁'),(-al*Q2,'−αQ₂ (α=0.8, opposite polarity)'),(S,'Q₁ − αQ₂ (net)')]):
 hp.mollview(m,sub=(2,3,i+1),title=t,min=-v,max=v,cmap='RdBu_r',cbar=False,hold=False)
ax=fig.add_subplot(2,1,2);th=np.linspace(0,np.pi,361)
for m,c,lab in[(Q1,'tab:red','Q₁'),(-al*Q2,'tab:blue','−αQ₂'),(S,'k','net')]:
 ax.plot(np.degrees(th),hp.get_interp_val(m,th,np.zeros_like(th)),color=c,lw=3 if c=='k' else 2,label=f'{lab}: rms={m.std():.3f}')
ax.set_xlabel('colatitude along φ=0 meridian [deg]');ax.set_ylabel('amplitude');ax.axhline(0,color='grey')
ax.legend(loc='upper center',bbox_to_anchor=(0.5,-0.2),ncol=3,frameon=False)
fig.suptitle('Illustration: two ℓ=2 patterns of opposite polarity superpose → reduced net quadrupole')
fig.tight_layout();fig.savefig('own_quadrupole_cancellation.png',dpi=120)
c1=lambda m:hp.anafast(m,lmax=2)[2];print(c1(Q1),c1(al*Q2),c1(S))
