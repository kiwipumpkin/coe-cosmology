import sys,os,json,numpy as np
sys.path.insert(0,'/workspace/clump_scan/sequential');os.chdir('/workspace/clump_scan/sequential')
import seq_core as S
os.chdir('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/ring_plane_cancel')
th=S.TH if hasattr(S,'TH') else np.radians(np.arange(360))
def mm(v):  # per-m rms in x CMB (1e-5); quadrature sum over m=1..12 + rest = total rms/FLOOR
    m=np.arange(1,61)[:,None];a=np.hypot(2*np.mean(v*np.cos(m*th),-1),2*np.mean(v*np.sin(m*th),-1));return a/np.sqrt(2)/S.FLOOR
f=np.load('/workspace/clump_scan/sequential_open/final_comp_s3.npz')
w0=S.pattern(f['X'],f['M']);w1=f['w'];R=S.R
out={}
for k,v in (('smooth_R',R),('web_start_residual',R-w0),('web_final_residual',R-w1),('web_final_pattern',w1)):
    a=mm(v);out[k]=dict(total=float(S.rms(v)/S.FLOOR),m=a[:12].tolist(),m13_60=float(np.sqrt(np.sum(a[12:]**2))));print(k,round(out[k]['total'],3),np.round(a[:12],3).tolist(),round(out[k]['m13_60'],3))
print('recomputed final left',S.rms(R-S.pattern(f['Xp'],f['Mp']))/S.FLOOR,'stored',float(f['left']))
np.savez('ring_patterns.npz',R=R,w0=w0,w1=w1,th=th);json.dump(out,open('decomp.json','w'),indent=1)
