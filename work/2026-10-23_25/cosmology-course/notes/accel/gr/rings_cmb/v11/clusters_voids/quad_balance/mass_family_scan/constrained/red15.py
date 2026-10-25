import os,json,numpy as np,warnings;warnings.filterwarnings('ignore')
H=os.getcwd();os.chdir('/workspace/cosmology-course/notes/accel/gr/three_vel/density')
exec(open('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/migkas_prof.py').read().split('rng=np.random')[0]);os.chdir(H)
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
D=json.load(open('l0scan_results.json'));e=[d for d in D if d['run']=='path10_RED' and d['l0']==15][0]
th=np.array(e['theta']);rat=e['ratio'];m=model(th,rat)
g=[k for k,v in list(globals().items()) if isinstance(v,np.ndarray) and v.shape==Y.shape and k not in('Y',)]
print(g)
print(L.min(),L.max(),S[:3],e['chi2'],e['R2'],th)
o=np.argsort(L%360);plt.figure(figsize=(11,6))
plt.errorbar(L%360,Y,yerr=S,fmt='ko',ms=3,label='Migkas 71 bins')
plt.plot((L%360)[o],m[o],'r-',lw=2,label='RED, l0=15 fixed')
plt.xlabel('galactic longitude (deg)');plt.ylabel('A/A_all');plt.xlim(0,360)
plt.title('RED universe, l0 = 15 deg fixed: chi2 = %.2f, R2 = %.3f'%(e['chi2'],e['R2']));plt.legend();plt.tight_layout();plt.savefig('red_migkas_l15.png',dpi=90)
