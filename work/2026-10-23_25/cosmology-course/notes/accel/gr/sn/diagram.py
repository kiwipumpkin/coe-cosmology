import json,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Wedge, FancyArrowPatch
RV=14.30048287260589
p=json.load(open('valley/valley_points.json'))[0]; v=[r for r in json.load(open('valley/valley_cmb.json')) if r['i']==0][0]
S=json.load(open('sn/sn_results.json')); Q={round(q['look_angle']):q for q in S['q_0']}
D,Rin,Rout=p['D'],p['Rin'],p['Rout']; rex=v['rmin']*RV; zstop=Q[179]['zmax_reached']
print('D',D,'Rin',Rin,'Rout',Rout,'excluded radius',rex,'core frac',p['mcc'],'zstop',zstop)
look={0:1,60:60,90:90,120:120,180:179}
def draw(ax,zoom):
    ax.add_patch(Wedge((0,0),Rout,0,360,width=Rout-Rin,color='#c9d9ee',label=f'uniform matter shell {Rin:.1f}–{Rout:.1f} Gpc'))
    ax.add_patch(Circle((0,0),rex,color='0.75',label=f'excluded inner region (radius {rex:.1f} Gpc)'))
    ax.plot(0,0,'k*',ms=16,label=f'CoE core ({p["mcc"]:.3f} of critical mass inside D)')
    ax.add_patch(Circle((D,0),RV,fill=False,ec='C2',lw=2.5,label=f'our visible sphere, radius {RV:.1f} Gpc'))
    ax.plot(D,0,'o',color='C3',ms=10,label=f'Earth, D = {D:.1f} Gpc')
    if not zoom: ax.text(0,(rex+Rin)/2,'empty gap (no matter)',ha='center',fontsize=12)
    for ang,c in zip(look,['C0','C4','C1','C4','C0']):
        a=np.radians(180-ang)  # ring angle 0 points away from CoE (+x)
        a=np.radians(ang); d=np.array([np.cos(a),np.sin(a)])
        if ang==180: d=np.array([-1,0])
        L=RV*0.95; ax.annotate('',xy=(D+L*d[0],L*d[1]),xytext=(D,0),arrowprops=dict(arrowstyle='->',color='C1',lw=1.8))
        q=Q[look[ang]]; txt=f"{ang}°: q0 {q['q0']:+.2f}, H0 ×{q['H0_over_global']:.3f}"+(f"\n(rays stop at z≈{zstop:.2f})" if ang==180 else "")
        if zoom: ax.text(D+1.08*L*d[0],1.08*L*d[1]+(3.0 if ang==0 else (1.2 if ang==180 else 0)),txt,fontsize=11,ha='center',va='bottom' if ang!=180 else 'bottom')
    s=1.0 if zoom else 1.0
    ax.annotate('',xy=(D+9,-RV-3),xytext=(D,-RV-3),arrowprops=dict(arrowstyle='-|>',color='k',lw=2.5))
    ax.text(D+4.5,-RV-5.5,'space carries us away from CoE\n(toy: 0.3c per ring radius)',ha='center',va='top',fontsize=11)
    ax.annotate('',xy=(D+1.2,1.4),xytext=(D,1.4),arrowprops=dict(arrowstyle='-|>',color='C3',lw=1.5))
    if zoom: ax.text(D+1.0,-2.6,'peculiar motion 0.00123c (arrow not to scale)',fontsize=10,color='C3')
    for a in (40,140,220,320):
        a=np.radians(a); ax.annotate('',xy=(D+(RV+3)*np.cos(a),(RV+3)*np.sin(a)),xytext=(D+RV*np.cos(a),RV*np.sin(a)),arrowprops=dict(arrowstyle='-|>',color='C2',lw=1.5))
    ax.text(D+RV*0.72+2,RV*0.72+2.2,'visible sphere edge\nrecedes at c',color='C2',fontsize=11)
    if zoom: ax.text(D+RV+1,-1.8,'0° (away from CoE,\ntoward Galactic l≈15°)',fontsize=11)
    if zoom: ax.text(D-RV-1,-3.0,'180°\n(toward CoE)',fontsize=11,ha='right')
    ax.set_aspect('equal'); ax.set_xlabel('distance along CoE axis (Gpc)'); ax.set_ylabel('Gpc')
plt.rcParams.update({'font.size':13})
fig,ax=plt.subplots(figsize=(14,10)); draw(ax,False); ax.set_xlim(-Rout-3,D+RV+16); ax.set_ylim(-Rout-3,Rout+3)
ax.text(D-RV,RV+3,'sight lines and labels: see zoom',fontsize=11)
ax.text(-Rout+1,-Rout+2,f"CMB ring RMS: model {v['rms']:.4f} vs observed 0.00087",fontsize=12,bbox=dict(fc='w',ec='k'))
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=2,fontsize=11,frameon=False); fig.tight_layout(); fig.savefig('scenario_diagram.png',dpi=150)
fig,ax=plt.subplots(figsize=(13,10)); draw(ax,True); ax.set_xlim(D-RV-12,D+RV+14); ax.set_ylim(-RV-10,RV+6)
ax.set_title('zoom on the visible sphere: supernova sight lines with inferred q0 and H0 (ratio to our shell)',fontsize=12)
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.04),ncol=2,fontsize=11,frameon=False); fig.tight_layout(); fig.savefig('scenario_diagram_zoom.png',dpi=150)
