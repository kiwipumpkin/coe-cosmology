import os;os.environ['GL']='1'
import sys,io,contextlib,numpy as np,matplotlib.pyplot as plt
src=open('pec.py').read().split("DL={0:")[0]   # compute part only
with contextlib.redirect_stdout(io.StringIO()): exec(src)
DL={0:'toward anti-CoE (0°, +15°)',90:'sideways toward (90°, 0°)',180:'toward CoE (180°, −15°)',270:'sideways toward (270°, 0°)'}
SEL=[(RUNS[2][1],'C2'),(RUNS[1][1],'C1')]
from scipy.signal import argrelextrema
def ext(y):
    yy=np.r_[y[-1],y,y[0]];mx=[i for i in range(1,len(yy)-1) if yy[i]>yy[i-1] and yy[i]>=yy[i+1]];mn=[i for i in range(1,len(yy)-1) if yy[i]<yy[i-1] and yy[i]<=yy[i+1]]
    return [round(float(XS[i-1]),1) for i in mx],[round(float(XS[i-1]),1) for i in mn]
def panel(ax,vd,fs=10):
    lo,hi=0,0
    for name,col in SEL:
        s=res[name][str(vd)];y=sel(curves[(name,vd)][1]);lo=min(lo,y.min());hi=max(hi,y.max())
        ax.plot(XS,y,color=col,lw=2.4,label=f"{name}: leftover {s['leftover_over_Planck']:.3f}× Planck (change {s['leftover_change_over_Planck']:+.3f})")
        pk,tr=ext(y);print(f"{DL[vd]:28s} {name:58s} peaks l={pk} troughs l={tr}")
    pad=0.12*(hi-lo);ax.set_ylim(lo-pad,hi+pad*2.2)
    planck(ax);ax.axvline({0:0.3,90:90,180:180,270:270}[vd],color='k',ls=':',lw=1)
    ax.set_xlim(0,360);ax.set_xticks(range(0,361,30));ax.grid(alpha=.25);ax.axhline(0,color='k',lw=.5)
    ax.set_ylabel('dT/T minus fixed no-motion axial monopole + cos dipole')
    ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),fontsize=fs,frameon=False,title=f'observer moving 370 km/s {DL[vd]} (dotted = motion direction)')
XLB='galactic longitude l (°) on great circle through anti-CoE (0°, +15°) and CoE (180°, −15°), tilted 15°; passes (90°, 0°), (270°, 0°)'
fig,axs=plt.subplots(2,2,figsize=(20,13))
for ax,vd in zip(axs.flat,DIRS): panel(ax,vd);ax.set_xlabel(XLB,fontsize=9)
fig.suptitle('Residual after removing only the fixed no-motion axial dipole (370 km/s motion dipole kept), speed 4c; Planck ring range %.7f'%PR,y=1.0,fontsize=13)
fig.tight_layout();fig.text(0.5,-0.03,CAPG,ha='center',fontsize=11);fig.savefig('pec_green_orange_2x2_gl.png',dpi=80,bbox_inches='tight');plt.close(fig)
fig,ax=plt.subplots(figsize=(16,9));panel(ax,270,12);ax.set_xlabel(XLB)
fig.text(0.5,-0.02,CAPG,ha='center',fontsize=10);fig.savefig('pec_green_orange_270_gl.png',dpi=90,bbox_inches='tight')
