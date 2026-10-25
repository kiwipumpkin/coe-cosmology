import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
import joint
j=json.load(open('joint.json'));x=j['x']
fig,(a1,a2)=plt.subplots(1,2,figsize=(20,9),gridspec_kw=dict(width_ratios=[1,1.3]))
a1.axis('off')
R2=1-j['chi2']/joint.CH0
rows=[['Constraint','Value','Pass?'],['Age ≥ 13 Gyr',f"{j['age']:.2f} Gyr",'~yes'],['CMB cancel ≈95.9%',f"{j['cancel']:.1f}% ({j['rms']/1e-5:.0f}× floor)",'~yes'],
['<H> in 67–73',f"{j['Hmean']:.1f}",'edge'],['Migkas χ² (null 45.1)',f"{j['chi2']:.0f}  (R² {R2:.1f})",'NO'],
['β₁, q, μ',f"{x[0]:.3f}, {x[1]:.3f}, {x[2]:.4f}",''],['s, φ, p',f"{x[3]:.2f}, {x[4]:.1f}°, {x[5]:.3f}",''],
['H⊥, D',f"{x[6]:.1f} km/s/Mpc, {j['D']/1e3:.2f} Gpc",''],['M',f"{j['M']:.2e} kg",''],['H(α) half-amp @350 Mpc',f"{100*j['Hamp']:.0f}%",'']]
t=a1.table(cellText=rows,loc='center',cellLoc='left');t.auto_set_font_size(False);t.set_fontsize(15);t.scale(1,3)
a1.set_title('Joint fit (v_o<c)',fontsize=22)
L=np.degrees(joint.LB);a2.errorbar(L,joint.Y,yerr=1/np.sqrt(joint.W),fmt='o',c='k',label='Migkas A/A_all')
a2.plot(L,j['A'],'r-',lw=3,label='Model H²/<H²> (350 Mpc)');a2.axhline(1,ls='--',c='gray',label='flat sky')
a2.set_xlabel('longitude (deg)',fontsize=18);a2.set_ylabel('A/A_all',fontsize=18);a2.tick_params(labelsize=15);a2.legend(fontsize=15)
plt.tight_layout();plt.savefig('joint_fit_summary.png',dpi=100)
