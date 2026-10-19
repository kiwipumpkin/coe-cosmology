import sys;sys.path.insert(0,'..');import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from free import resid,vec,lb,PR
from scipy.spatial.transform import Rotation as R
F=json.load(open('free_fit.json'));S0=json.load(open('../rotfit_stats.json'))
e1,e2=vec(0,15),vec(90,0);c=np.radians(np.arange(0,360,.5));N=np.outer(np.cos(c),e1)+np.outer(np.sin(c),e2)
L=np.degrees(np.arctan2(N[:,1],N[:,0]))%360;o=np.argsort(L);L=L[o];N=N[o]
a0=vec(0,15);v0=vec(*S0['fit']['motion_lb'])
def curve(a,v):
    import fit;Vs=fit.V;fit.V=N
    import free;free.V=N
    from fit import sky
    m=sky(a,v);r=m-free.c0-free.c1*(N@a);fit.V=Vs;free.V=Vs;return r
def pk(y):return L[np.argmax(y)],L[np.argmin(y)]
def dl(x,y):return (x-y+180)%360-180
groups=[('yaw: rotate about the galactic pole (shift in l)','z',[0,45,90,180]),('tilt about the axis through l = 0° (in the galactic plane)','x',[0,15,30,45]),('tilt about the axis through l = 90° (in the galactic plane)','y',[0,15,30,45])]
fig,axs=plt.subplots(2,2,figsize=(22,15));rep={}
for ax,(t,axn,steps) in zip(axs.flat,groups):
    h0,c0_=None,None
    for s,col in zip(steps,['C1','C0','C2','C3']):
        Rm=R.from_euler(axn,s,degrees=True).as_matrix();y=curve(Rm@a0,Rm@v0);h,cl=pk(y)
        if s==0:h0,c0_=h,cl
        rep[f'{axn}{s}']=dict(hot_l=float(h),cold_l=float(cl),hot_shift=float(dl(h,h0)),cold_shift=float(dl(cl,c0_)))
        ax.plot(L,y*1e3,col,lw=2.4,label=f"{s}°: hot peak l = {h:.1f}° (moved {dl(h,h0):+.1f}°), cold l = {cl:.1f}° (moved {dl(cl,c0_):+.1f}°)")
    ax.set_title(t,fontsize=15,pad=118)
for ax,(lab,x,col) in zip([axs[1,1]]*3,[('paper axis CoE (180°, −15°)',[0,15,*S0['fit']['motion_lb']],'C1'),('best fit A, CoE (264.0°, +48.2°)',F['A_full_dipole']['x'],'C0'),('best fit B, CoE (77.3°, −57.9°)',F['B_axial_dipole_removed']['x'],'C3')]):
    y=curve(vec(x[0],x[1]),vec(x[2],x[3]));h,cl=pk(y);rep[lab]=dict(hot_l=float(h),cold_l=float(cl))
    ax.plot(L,y*1e3,col,lw=2.4,label=f"{lab}: hot l = {h:.1f}°, cold l = {cl:.1f}°")
axs[1,1].set_title('fitted axes (each with its own fitted motion)',fontsize=15,pad=118)
for ax in axs.flat:
    ax.axvline(264.02,color='r',lw=2);ax.axvline(84.02,color='b',lw=2);ax.text(264.02,ax.get_ylim()[1],' Planck hot 264.02°',color='r',va='top',fontsize=12);ax.text(84.02,ax.get_ylim()[1],' Planck cold 84.02°',color='b',va='top',fontsize=12)
    ax.set_xlim(0,360);ax.set_xticks(range(0,361,30));ax.grid(alpha=.3);ax.tick_params(labelsize=12)
    ax.set_ylabel('dT/T minus axial dipole (×10⁻³)',fontsize=13);ax.set_xlabel('galactic longitude l (°) on the great circle through (0°, +15°), (90°, 0°), (180°, −15°), (270°, 0°)',fontsize=11)
    ax.legend(loc='lower center',bbox_to_anchor=(.5,1.0),fontsize=12,frameon=False)
fig.suptitle('Own crossing, distributed mass b, 370 km/s motion toward (254.7°, +46.2°) rotated with the sky; peaks = highest/lowest point on this circle',fontsize=15,y=1.0)
fig.tight_layout();fig.savefig('rotation_effect_2x2.png',dpi=75,bbox_inches='tight');json.dump(rep,open('peak_shifts.json','w'),indent=1);print(json.dumps(rep,indent=0))
