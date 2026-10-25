exec(open('sp.py').read().split('h0,c0=pk')[0])
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
h0,c0=pk(curve(a0,v0));out={}
def draw(ax,items,title):
    for lab,Rm,col in items:
        y=curve(Rm@a0,Rm@v0);h,cl=pk(y);ax.plot(L,y*1e3,col,lw=2.4,label=f"{lab}: hot l = {h:.1f}° ({dl(h,h0):+.1f}°), cold l = {cl:.1f}° ({dl(cl,c0):+.1f}°)");out[title+' | '+lab]=dict(hot=float(h),cold=float(cl),dhot=float(dl(h,h0)),dcold=float(dl(cl,c0)))
    ax.set_title(title,fontsize=15,pad=100)
C=['C0','C1','C2','C3']
fig,axs=plt.subplots(2,2,figsize=(22,15))
draw(axs[0,0],[(f'{t:+d}°',rot('y',t),c) for t,c in zip([-10,-5,5,10],C)],'tilt about the axis through l = 90°, −10° to +10° (no tilt: hot 255.5°, cold 174.6°)')
draw(axs[0,1],[(f'{t:+d}°',rot('y',t),c) for t,c in zip([0,30,60,90],C)],'same tilt, larger angles (where the spread becomes real)')
draw(axs[1,0],[(f'{t:+d}°',rot('x',t),c) for t,c in zip([-10,-5,5,10],C)],'for contrast: tilt about the axis through l = 0°, −10° to +10° (pinch)')
Rb=rot('zxy',[15,45,90])
draw(axs[1,1],[('paper axis (no rotation)',np.eye(3),'C1'),('yaw 15°, tilt 45° about l = 0°, tilt 90° about l = 90°',Rb,'C3')],'best combination found (grid of 15° steps)')
for ax in axs.flat:
    y1=ax.get_ylim()[1];ax.axvline(264.02,color='r',lw=2);ax.axvline(84.02,color='b',lw=2);ax.text(264.02,y1,' Planck hot 264.02°',color='r',va='top',fontsize=12);ax.text(84.02,y1,' Planck cold 84.02°',color='b',va='top',fontsize=12)
    ax.set_xlim(0,360);ax.set_xticks(range(0,361,30));ax.grid(alpha=.3);ax.tick_params(labelsize=12);ax.set_ylabel('dT/T minus axial dipole (×10⁻³)',fontsize=13)
    ax.set_xlabel('galactic longitude l (°) on the great circle through (0°, +15°), (90°, 0°), (180°, −15°), (270°, 0°)',fontsize=11);ax.legend(loc='lower center',bbox_to_anchor=(.5,1.0),fontsize=12,frameon=False)
fig.suptitle('Own crossing, distributed mass b, 370 km/s motion (254.7°, +46.2°) rotated with the sky; shifts in brackets are from the unrotated curve',fontsize=15,y=1.0)
fig.tight_layout();fig.savefig('spread_2x2.png',dpi=75,bbox_inches='tight')
# full-sky check of best combination
import healpy as hp
a=Rb@a0;v=Rb@v0;m=fit.sky(a,v);d=fit.dip(m);q=fit.quad(m)[0];r=free.resid(a,v);dr=fit.dip(r);qr=fit.quad(r)[0]
PH,PQ=fit.PH,fit.PQ;h=lambda x,y:min(fit.ang(x,y),180-fit.ang(x,y))
out['best']=dict(antiCoE=lb(a),CoE=lb(-a),motion=lb(v),full=dict(dip_hot=lb(d),hot_off=float(fit.ang(d,PH)),cold_off=float(fit.ang(-d,-PH)),quad=lb(q),quad_off=float(h(q,PQ)),amp_over_Planck=float(np.linalg.norm(d)*fit.T0/3362)),
  residual=dict(dip_hot=lb(dr),hot_off=float(fit.ang(dr,PH)),quad=lb(qr),quad_off=float(h(qr,PQ)),amp_over_Planck=float(np.linalg.norm(dr)*fit.T0/3362)))
json.dump(out,open('spread_plot.json','w'),indent=1,default=float);print(json.dumps(out['best'],default=float))
