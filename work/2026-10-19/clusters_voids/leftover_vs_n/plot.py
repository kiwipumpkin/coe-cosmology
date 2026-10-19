import json,re,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
def logc(f):
    c=[(0,5.435961527713568)];[c.append((int(m[0]),float(m[1]))) for m in re.findall(r'rms (\d+) ([\d.]+)',open(f).read())];return np.array(c)
D={f:np.array(json.load(open(f'res_{f}_dip.json'))['curve']) for f in ('mass','void','both')}
R={f:logc(f'log_{f}_rms.txt') for f in ('mass','void','both')}
ext={}
for f,c in R.items():
    k,b=np.polyfit(c[-8:,0],c[-8:,1],1);ext[f]=float((1-b)/k) if k<0 else None
fig,axs=plt.subplots(1,2,figsize=(22,9))
a=axs[0]
for f,col,lab in (('mass','#cb181d','masses only'),('void','#3182bd','voids only'),('both','#6a3d9a','masses + voids')):
    c=D[f];a.plot(c[:,0],c[:,1],'o-',color=col,lw=3,ms=10,label=f"on the recombination surface, {lab}: 1× at N = {int(c[c[:,1]<=1][0,0])}")
a.plot([0,6],[1.722,1.722],'--',color='.4',lw=3,label='inside or outside our sphere, not touching the surface: no change (1.72×)')
a.axhline(1,color='k',lw=2);a.text(4.6,1.05,'1× target',fontsize=14)
a.set_title('(a) dipole only (× Planck dipole)',fontsize=17);a.set_xlabel('number of bodies N',fontsize=16);a.set_ylabel('leftover (× Planck dipole)',fontsize=16)
a=axs[1]
for f,col,lab in (('mass','#cb181d','masses only'),('void','#3182bd','voids only'),('both','#6a3d9a','masses + voids')):
    c=R[f];a.plot(c[:,0],c[:,1],'o-',color=col,lw=3,ms=8,label=f"surface, {lab}: {c[-1,1]:.2f}× at N = {int(c[-1,0])} (stopped); straight-line trend reaches 1× near N ≈ {ext[f]:.0f}")
a.plot([0,25],[5.436,5.436],'--',color='.4',lw=3,label='inside or outside, not touching the surface: no change (5.44×)')
a.axhline(1,color='k',lw=2);a.set_ylim(0,6.5)
a.set_title('(b) all patterns l = 1–12 (RMS, × Planck dipole RMS)',fontsize=17);a.set_xlabel('number of bodies N',fontsize=16)
for a in axs:a.grid(alpha=.3);a.tick_params(labelsize=14);a.legend(loc='lower center',bbox_to_anchor=(.5,1.06),fontsize=12,frameon=False)
fig.text(.5,-.05,'Our model: smooth shrinking-radius sky (mass b, CatWISE slope, 1+z ≈ 1100). A body changes the CMB only where it cuts the recombination surface (taken at 14.3 Gpc): dT/T = ln(1+δ)/3 over that patch.\nBody sizes and contrasts limited to Chris\'s seed-3 laws (clusters ≤ 0.30 Gpc, contrast 30–300; voids ≤ 2.14 Gpc, depth −0.45 to −0.25). Greedy addition with joint refit after each body.',ha='center',fontsize=12)
fig.savefig('leftover_vs_N.png',dpi=80,bbox_inches='tight')
json.dump(dict(dipole={f:D[f].tolist() for f in D},rms={f:R[f].tolist() for f in R},rms_linear_trend_N_for_1x=ext,rms_after_dipole_cancel={'mass_4':12.826,'void_1':12.899}),open('leftover_vs_N.json','w'),indent=1)
print(ext,{f:R[f][-1] for f in R})
