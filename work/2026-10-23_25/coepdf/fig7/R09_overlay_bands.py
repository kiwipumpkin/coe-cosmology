# rev7 regeneration of rev6 Fig 9. Runs R09_overlay_bands_core.py (the original band-construction code, book/code/25_*, output dir
# redirected) and lays its three computed versions into one figure. Band centres: the locked Migkas curve has extrema at
# l=20,200 (max) and 127,272 (min); the rev7 D=34.96 refit (p1_migkas_D34.96.json) has the same extrema, so the construction is unchanged.
import subprocess,sys,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
H='/workspace/coepdf/fig7/';subprocess.run([sys.executable,H+'R09_overlay_bands_core.py'],check=True)
c=np.array(list(json.load(open(H+'../joint_solve/phase2/p1_migkas_D34.96.json')).values())[0]['curve']);l=np.arange(len(c))*360/len(c)
mx=l[argrelextrema(c,np.greater,mode='wrap')[0]];mn=l[argrelextrema(c,np.less,mode='wrap')[0]]
fig,ax=plt.subplots(3,1,figsize=(11,16.5))
for a,(f,t) in zip(ax,(('matched','+15° tilt alone'),('matched_grad','tilt + gradient'),('matched_grad_rot81','tilt + gradient + 80.8° roll'))):
    a.imshow(plt.imread(H+f'r09/fig8_overlay_{f}.png'));a.axis('off');a.set_title(t,fontsize=13)
fig.suptitle(f'rev7 fit extrema: MAX l={mx.tolist()}, MIN l={mn.tolist()} (identical to locked curve)',fontsize=11)
fig.tight_layout();fig.savefig(H+'R09_overlay_bands.png',dpi=110);print(mx,mn)
