import os,sys,json,numpy as np
H=os.getcwd()
import interp_an as I
import matplotlib.pyplot as plt
R=json.load(open(H+'/interp_results.json'));PL=np.loadtxt('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb/planckTT.txt')
def fig(k,name):
    A=R[k]
    if 'chi2' not in A:print('skip',k,A);return
    f,ax=plt.subplots(1,3,figsize=(19,5.2))
    a=ax[0];a.errorbar(I.L,I.Y,yerr=I.S,fmt='o',ms=3,color='k',label='Migkas 71 bins');a.plot(I.Ls,I.curve(np.array(A['theta']),A['ratio']),'C3',lw=1.8,label='model')
    a.set_xlim(0,360);a.set_xlabel('galactic longitude [deg]');a.set_ylabel('A/A_all');a.legend(fontsize=7)
    a.set_title('χ²=%.2f  R²=%.3f  ℓ0=%.1f°%s  ratio=%.5f'%(A['chi2'],A['R2'],A['l0'],' (at bound)' if A['l0_at_bound'] else '',A['ratio']),fontsize=9)
    o=np.load(A['run'])['o'];o=o[np.argsort(o[:,0])];psi=o[:,0];lz=np.log(o[:,1]);w=np.sin(np.radians(psi))
    d=-(lz-np.trapezoid(lz*w,psi)/np.trapezoid(w,psi));mu=np.cos(np.radians(psi));nd=d-A['a'][0]*mu
    a=ax[1];a.plot(psi,d,'C3o-',ms=2,label='monopole removed');a.plot(psi,nd,'C0s-',ms=2,label='monopole + dipole (a1 P1) removed')
    a.invert_xaxis();a.set_xlabel('ψ [deg] (180 = toward CoE, left)');a.set_ylabel('ΔT/T');a.legend(fontsize=7)
    a.set_title('a1..a4 = '+', '.join('%.2e'%x for x in A['a']),fontsize=9)
    a=ax[2];ll=np.arange(2,31);a.semilogy(ll,np.maximum(np.array(A['Dl'][2:31]),1e-8),'C3o-',ms=3,label='model')
    a.errorbar(PL[:29,0],PL[:29,1],yerr=[PL[:29,2],PL[:29,3]],fmt='ks',ms=3,capsize=2,label='Planck 2018 TT');a.set_xlabel('ℓ');a.set_ylabel('D_ℓ [μK²]');a.legend(fontsize=7)
    a.set_title('D2..D6 = '+', '.join('%.3g'%x for x in A['Dl'][2:7]),fontsize=9)
    f.suptitle('%s   C=%s  frac=%g'%(k,A['C'],A['frac']),fontsize=10);plt.tight_layout();plt.savefig(H+'/'+name,dpi=75);plt.close()
for k in range(11):
    key=[x for x in R if x.startswith('path%02d'%k)][0];fig(key,'interp_step%02d_panels.png'%k)
for k in [x for x in R if x.startswith('hyb')]:fig(k,'%s_panels.png'%k)
st=[];D2=[];D3=[];C2=[]
for k in range(11):
    A=R[[x for x in R if x.startswith('path%02d'%k)][0]]
    if 'chi2' in A:st.append(k);D2.append(A['Dl'][2]);D3.append(A['Dl'][3]);C2.append(A['chi2'])
f,ax=plt.subplots(1,2,figsize=(12,4.5));ax[0].semilogy(st,D2,'o-',label='D2');ax[0].axhline(226,color='r',ls='--',label='Planck D2=226 μK²');ax[0].set_xlabel('path step');ax[0].set_ylabel('D2 [μK²]');ax[0].semilogy(st,D3,'s-',label='D3');ax[0].axhline(937,color='m',ls=':',label='Planck D3=937 μK²');ax[0].set_ylabel('D_ℓ [μK²]');ax[0].legend()
ax[1].plot(st,C2,'o-');ax[1].set_xlabel('path step');ax[1].set_ylabel('Migkas χ²');plt.tight_layout();plt.savefig(H+'/interp_summary.png',dpi=80)
