# Axial-mode analysis: rotate a_lm so the CoE axis (Galactic l=180, b=-15) is the z-axis; m=0 power fraction vs l.
import json,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.signal import find_peaks
LM=60;T0=2.7255e6;th0=np.radians(90+15);ph0=np.radians(180)
def rot(alm):a=alm.copy();hp.rotate_alm(a,-ph0,-th0,0.0,lmax=LM);return a   # zyz: bring (th0,ph0) to the pole
# self-test: axisymmetric blob at the CoE must be pure m=0 after rotation
v=hp.ang2vec(th0,ph0);ns=64;Vp=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T;tm=np.exp(-np.arccos(np.clip(Vp@v,-1,1))**2/0.1)
def f0(m):
    a=rot(hp.map2alm(m,lmax=LM));l=np.arange(LM+1);out=[];P0=[]
    for L in l:
        i0=hp.Alm.getidx(LM,L,0);tot=abs(a[i0])**2+2*sum(abs(a[hp.Alm.getidx(LM,L,mm)])**2 for mm in range(1,L+1))
        out.append(abs(a[i0])**2/tot if tot>0 else 0);P0.append(L*(L+1)/(2*np.pi)*abs(a[i0])**2/(2*L+1))
    return np.array(out),np.array(P0)
print('self-test m0 fraction l=1..6:',np.round(f0(tm)[0][1:7],4))
z=np.load('/workspace/redweb2/realweb/los_flow_R/hl_R4_ns64.npz');bg=z['sm']+z['rs_cc'];fld=np.load('mc_sky0.npy')
sm=hp.ud_grade(hp.read_map('/workspace/planck/smica.fits'),256)/1e0
S={'Planck SMICA (unmasked, K)':hp.remove_dipole(sm),'CoE background (smooth red+RS_cc)':hp.remove_dipole(bg),'stochastic field, 1 realization':fld,'field+background':hp.remove_dipole(fld+bg)}
l=np.arange(LM+1);R={};fig,ax=plt.subplots(2,1,figsize=(11,9))
for j,(k,m) in enumerate(S.items()):
    f,P0=f0(m);exc=f*(2*l+1);R[k]=dict(f0=f.tolist(),excess=exc.tolist())
    pk,_=find_peaks(exc[2:31],height=3);pk=pk+2;R[k]['peaks_l2_30_excess>3']=pk.tolist()
    lad=[(a,b,c) for a in pk for b in pk for c in pk if a<b<c and abs(b/a-2)<0.15 and abs(c/a-3)<0.15]
    R[k]['ladders_1_2_3']=[list(map(int,t)) for t in lad]
    R[k]['implied_L_Gpc_R4']=[float(2*np.pi*4.0/t[0]) for t in lad];R[k]['implied_L_Gpc_R4.23']=[float(2*np.pi*4.23/t[0]) for t in lad]
    ax[0].plot(l[2:],exc[2:],'-o',ms=3,label=k)
    print(k,'| m0 excess l2-10:',np.round(exc[2:11],2),'| peaks>3:',pk.tolist(),'| 1:2:3 ladders:',lad)
ax[0].axhline(1,color='k',ls=':',label='isotropic expectation');ax[0].set_yscale('log');ax[0].set_ylabel('(2ℓ+1)|a_ℓ0|²/Σ|a_ℓm|²  (CoE frame)');ax[0].legend(loc='lower center',fontsize=7,ncol=2);ax[0].set_title('axial (m=0) power excess, CoE axis = z')
# ensemble expectation for isotropic Gaussian: f0*(2l+1) ~ chi2_1 / mean -> 95% level ~ 3.84
ax[0].axhline(3.84,color='gray',ls='--',label='95% (isotropic, per ℓ)')
fE=[np.array(R[k]['excess']) for k in S];ax[1].plot(l[2:],fE[0][2:],'k-o',ms=3,label='Planck');ax[1].plot(l[2:],fE[3][2:],'C3-o',ms=3,label='model field+bg')
for lj in (1,2,3):pass
ax[1].axhline(3.84,color='gray',ls='--');ax[1].set_xlabel('ℓ');ax[1].set_ylabel('m=0 excess');ax[1].legend(loc='lower center');ax[1].set_yscale('log')
fig.tight_layout();fig.savefig('f_axial_m0.png',dpi=100);json.dump(R,open('axial.json','w'),indent=1)
