# FIXED 3:3x PM: curve = squared-refit H/<H> (migkas_refit_sq.py), old CSV was H-like so sqrt was double-applied.
# NEW phase3 (H0 from sourced table h0_sources.json): mean=inverse-variance 68.64, low/high = band edges 66.90/74.08.
# NEW rev7: directional H0 implied by the locked Migkas A/A_all curve, with the sky mean fixed to the observed H0.
# Migkas et al. (2020, A&A 636, A15): L_X = 4 pi D_L^2 F and D_L ~ 1/H0, with H0 assumed isotropic the fitted L_X is low where true H0 is low (D_L too small), so A_apparent ~ H0^2 -> H0(l)/<H0> = (A/A_all)^(+1/2) (Migkas+2021 relation).
# The linear form H0 = <H0> * A/A_all (power wrong) is also reported for comparison only.
import numpy as np,json,matplotlib;J=json.load(open('/workspace/coepdf/joint_solve/phase3/h0_sources.json'));HM,HL,HH=round(J['H0_mean_invvar'],2),round(J['H0_low'],2),round(J['H0_high'],2);matplotlib.use('Agg');import matplotlib.pyplot as plt
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,y=d['l_deg'],d['A_over_Aall'];el,eh=d['A_err_low'],d['A_err_high']
c=np.loadtxt('/workspace/coepdf/joint_solve/phase3/migkas_refit_sq_curve.csv',delimiter=',',skiprows=1);lc,A,Hc=c[:,0],c[:,1],c[:,2]  # squared refit: A=af^2/<af^2>, Hc=af/<af> (H-like, no extra sqrt)
out={}
fig,ax=plt.subplots(1,2,figsize=(13,4.6),gridspec_kw=dict(width_ratios=[2.2,1]))
for H,ls,lab in ((HM,'-',f'inv-var mean {HM:g}'),(HL,':',f'band low {HL:g}'),(HH,'--',f'band high {HH:g}')):
    h=Hc.copy();h=h/h.mean()*H                      # sky-mean (longitude-mean) constraint
    hd=y**0.5;hd=hd/np.mean(hd)*H
    a1=2*np.mean(h*np.exp(-1j*np.radians(lc)));dip=abs(a1)/H;ldip=np.degrees(np.angle(a1))%360
    lin=Hc/Hc.mean()*H
    out[str(H)]=dict(mean=float(h.mean()),rms_frac=float(h.std()/H),min=float(h.min()),l_min=float(lc[h.argmin()]),max=float(h.max()),l_max=float(lc[h.argmax()]),
        peak_to_peak_frac=float((h.max()-h.min())/H),dipole_frac=float(dip),l_dipole_max=float(ldip),
        data_points=dict(min=float(hd.min()),max=float(hd.max()),rms_frac=float(hd.std()/H)),
        linear_form=dict(min=float(lin.min()),max=float(lin.max()),rms_frac=float(lin.std()/H)))
    ax[0].plot(lc,h,ls,color='C3' if H==HM else '0.4',lw=2 if H==HM else 1.2,label=f'squared-refit CoE curve, <H0>={H:g} ({lab})')
    if H==HM:
        ax[0].errorbar(l,hd,[hd-np.mean(y**0.5)**-1*H*(y-el)**0.5,np.mean(y**0.5)**-1*H*(y+eh)**0.5-hd],fmt='ko',ms=3,capsize=1.5,label='Migkas points, same transform')
        ax[1].bar(['rms','dipole','max-min'],[100*out[str(HM)]['rms_frac'],100*dip,100*out[str(HM)]['peak_to_peak_frac']],color='C3',label='CoE squared refit')
ax[0].axhspan(HL,HH,color='gold',alpha=.25,label=f'observed band {HL:g}-{HH:g} (span of Planck, DESI, CCHP, SH0ES ±1σ)')
ax[0].axvline(280,color='b',lw=.8);ax[0].text(282,ax[0].get_ylim()[0]+0.3,'Migkas 2021\nlow-H0 region l~280°',fontsize=7,color='b')
ax[0].set_xlim(0,360);ax[0].set_xlabel('Galactic longitude [deg]');ax[0].set_ylabel('H0(l) [km/s/Mpc]');ax[0].legend(fontsize=7.5,loc='upper right');ax[0].grid(alpha=.3)
ax[0].set_title('Directional H0 = <H0> (A/A_all)^(1/2), sky mean fixed',fontsize=10)
ax[1].axhline(9,color='b',ls='--',label='Migkas+2021: ~9% H0 variation')
ax[1].axhline(16,color='g',ls=':',label='Migkas+2020: max D_L variation 16±3%')
ax[1].set_ylabel('% of <H0>');ax[1].legend(fontsize=7.5);ax[1].grid(alpha=.3,axis='y');ax[1].set_title('Size of the anisotropy',fontsize=10)
fig.tight_layout();fig.savefig(__file__.replace('.py','.png'),dpi=140,facecolor='w');json.dump(out,open(__file__.replace('.py','.json'),'w'),indent=1)
print(json.dumps(out[str(HM)],indent=1))
