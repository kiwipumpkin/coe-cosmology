# Percent of l=2..30 power per l: light-path model (monopole+dipole removed map) vs Planck TT.
import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
m=np.load('cmb_lightpath_maps.npz');p=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt')
l=np.arange(2,31);c=m['cl_nodip'][2:31]*(2*l+1);pm=100*c/c.sum()
Lp=p[:,0].astype(int);s=(Lp>=2)&(Lp<=30);P=(2*Lp[s]+1)*2*np.pi*p[s,1]/(Lp[s]*(Lp[s]+1));pp=100*P/P.sum()
fig,ax=plt.subplots(figsize=(10,4.5));ax.bar(l-.2,pm,.4,label='model',color='C3');ax.bar(l+.2,pp,.4,label='Planck TT',color='k')
ax.set_xlabel('$\\ell$');ax.set_ylabel('% of $\\sum_{2}^{30}(2\\ell+1)C_\\ell$');ax.legend();ax.set_title('Share of l=2-30 power')
plt.savefig('cmb_lightpath_pct_l2_30.png',dpi=120,bbox_inches='tight')
json.dump(dict(l=l.tolist(),model_pct=pm.tolist(),planck_pct=pp.tolist()),open('cmb_lightpath_pct_l2_30.json','w'),indent=0)
print('model l2,3,4:',pm[:3],' planck:',pp[:3])
