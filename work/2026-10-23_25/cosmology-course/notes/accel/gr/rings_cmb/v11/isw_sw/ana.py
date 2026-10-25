import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from numpy.polynomial.legendre import legval
from scipy.interpolate import CubicSpline
x=np.linspace(-1,1,8001);ps=np.degrees(np.arccos(x));TC=2.7255e6
def a2(f):
    c=[0,0,1];return 5/2*np.trapezoid(f*legval(x,c),x)
res={};fig,ax=plt.subplots(1,2,figsize=(16,6.5))
for f,lab in [('r30.npz','M/Mcrit=30 (soft p=2)'),('r4.npz','M/Mcrit=4.362 (soft p=2)')]:
    d=np.load(f);o=d['o'];sp=lambda col:CubicSpline(o[:,0],col)(ps)
    geo=-sp(np.log(o[:,1]));geo-=geo.mean()
    sw=sp(o[:,4])/3;sw-=sw.mean()
    isw=sp(o[:,5]);isw-=isw.mean()
    P2=legval(x,[0,0,1]);g2,s2,i2=a2(geo),a2(sw),a2(isw)
    res[lab]=dict(geometric_a2=g2,SW_a2=s2,ISW_a2=i2,SWplusISW_a2=s2+i2,
      geo_axis_uK=g2*TC,SW_axis_uK=s2*TC,ISW_axis_uK=i2*TC,
      SW_over_geo=s2/g2,ISW_over_geo=i2/g2,total_over_geo=(g2+s2+i2)/g2,
      Phi_emit_range=[float(o[:,4].min()),float(o[:,4].max())],Phi_obs=float(d['Phi_obs']),
      ISW_range=[float(o[:,5].min()),float(o[:,5].max())],max_abs_Phi_over_c2=float(abs(o[:,4]).max()))
    k=0 if '30' in f else 1
    for v,c,n in [(g2,'k','geometric (ray-traced)'),(s2,'tab:blue','SW Φ_e/3'),(i2,'tab:red','ISW 2∫∂tΦ dt')]:
        ax[k].plot(ps,v*P2*TC,c=c,lw=2,label=f'{n}: a2={v:.3g}')
    ax[k].set_yscale('symlog',linthresh=1e4);ax[k].set_xlabel('ψ from CoE axis [deg]  (ψ=0 toward CoE, left)');ax[k].set_ylabel('l=2 component ΔT [μK]')
    ax[k].invert_xaxis();ax[k].axhline(0,c='gray',lw=.5);ax[k].set_title(lab);ax[k].legend(fontsize=9)
fig.suptitle('Naive linear Newtonian SW/ISW quadrupole vs geometric l=2 — NOTE |Φ|/c² ~1e5–1e6 at emission: linear SW/ISW invalid, and LTB geodesic already contains these effects')
plt.savefig('isw_sw_quad_estimate.png',dpi=110,bbox_inches='tight')
json.dump(res,open('isw_sw_quad_estimate.json','w'),indent=1);print(json.dumps(res,indent=1))
