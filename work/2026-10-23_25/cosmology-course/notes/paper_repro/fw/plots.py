import json,os,numpy as np,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
import core
L=np.arange(360.)
def mig(ax,fw,col,lab):
    m=json.load(open(f'migkas_{fw}.json'))['1deg']; c=np.array(m['curve'])
    ax.plot(L,c,col,lw=2.5,label=f"{lab}: R²={m['R2']:.3f}, χ²={m['chi2']:.2f}, {m['sigma']:.1f}σ")
    return m
def migdata(ax):
    ax.errorbar(core.LDEG,core.Y,yerr=[core.MLO,core.MHI],fmt='o',color='k',ms=4,alpha=.6,label='Migkas cluster data')
    ax.set_xlim(0,360); ax.set_xlabel('Galactic longitude (degrees)'); ax.set_ylabel('A / A_all'); ax.grid(alpha=.3)
def cmbcurves(fw):
    if not os.path.exists(f'cmb_{fw}_curves.json'): return {}
    C=json.load(open(f'cmb_{fw}_curves.json')); p=json.load(open(f'migkas_{fw}.json'))['1deg']['params']
    out={}
    for k,T in C.items():
        T=np.array(T); d=(T-T.mean())/T.mean()/1e-5; l=(np.arange(360)+p['l0'])%360; o=np.argsort(l); out[k]=(l[o],d[o])
    return out
LBL={'with_mass_f1.0_without_gravity':'With mass, all motion is space growing (f=1). Gravity left out: its factor is negative',
     'with_mass_f0.9_stretch_only':'With mass, f=0.9: stretch only (motion-through-space and gravity factors negative)',
     'no_mass_f1.0_all_effects':'No mass, f=1: all effects (motion-through-space = 0, gravity = 0)',
     'no_mass_f0.9_stretch_only':'No mass, f=0.9: stretch only (motion-through-space factor negative)'}
STY={'with_mass_f1.0_without_gravity':'r-','with_mass_f0.9_stretch_only':'r--','no_mass_f1.0_all_effects':'b-','no_mass_f0.9_stretch_only':'b--'}
YL=(-14000,14000)
# per framework
for fw,name in (('classical','Classical (Newton)'),('sr','Special relativity'),('gr','General relativity (LTB)')):
    if not os.path.exists(f'migkas_{fw}.json'): continue
    fig,(a1,a2)=plt.subplots(2,1,figsize=(12,10))
    migdata(a1); mig(a1,fw,'g',name+' fit'); a1.legend(fontsize=9); a1.set_title(name+': refit to Migkas, and CMB temperature vs direction')
    cc=cmbcurves(fw); n=0
    for k,(l,d) in cc.items():
        if k in LBL: a2.plot(l,d,STY[k],lw=2,label=LBL[k]); n+=1
    if fw=='gr' and cc:
        for k,(l,d) in cc.items(): a2.plot(l,d,'r-' if k=='with_mass' else 'b-',lw=2,label=('With mass' if k=='with_mass' else 'No mass')+': full light-path redshift'); n+=1
    if n==0: a2.text(180,0,{'sr':'No CMB can be computed: no light from the recombination radius reaches us now.\nMatter moves at 0.0024 of light speed in this fit, so that light crossed\nthe whole universe long ago (see report).',
                            'gr':'GR CMB not finished (see report).','classical':''}[fw],ha='center',va='center',fontsize=11)
    a2.axhspan(-1,1,color='g',alpha=.3,label='Required flatness: 1 part in 100,000')
    a2.set_xlim(0,360); a2.set_ylim(*YL); a2.set_xlabel('Galactic longitude (degrees)'); a2.set_ylabel('Temperature change, parts per 100,000'); a2.grid(alpha=.3); a2.legend(fontsize=8,loc='lower left')
    fig.tight_layout(); fig.savefig(f'cmb_{fw}.png',dpi=130); plt.close(fig)
# comparison
fig,(a1,a2)=plt.subplots(2,1,figsize=(12,10)); migdata(a1)
for fw,col,name in (('classical','g','Classical'),('sr','m--','Special relativity'),('gr','c:','General relativity')):
    if os.path.exists(f'migkas_{fw}.json'): mig(a1,fw,col,name)
a1.legend(fontsize=9); a1.set_title('All three frameworks: Migkas fit (top) and CMB vs direction (bottom), same scales')
for fw,tag in (('classical','Classical'),('gr','GR')):
    for k,(l,d) in cmbcurves(fw).items():
        if fw=='classical' and k in ('with_mass_f1.0_without_gravity','no_mass_f1.0_all_effects'): a2.plot(l,d,STY[k],lw=2,label='Classical: '+LBL[k])
        if fw=='gr': a2.plot(l,d,'m-' if k=='with_mass' else 'c-',lw=2,label='GR: '+k.replace('_',' '))
a2.axhspan(-1,1,color='g',alpha=.3,label='Required flatness: 1 part in 100,000'); a2.text(5,11000,'Special relativity: no CMB computable',fontsize=10)
a2.set_xlim(0,360); a2.set_ylim(*YL); a2.set_xlabel('Galactic longitude (degrees)'); a2.set_ylabel('Temperature change, parts per 100,000'); a2.grid(alpha=.3); a2.legend(fontsize=8,loc='lower left')
fig.tight_layout(); fig.savefig('cmb_compare.png',dpi=130)
