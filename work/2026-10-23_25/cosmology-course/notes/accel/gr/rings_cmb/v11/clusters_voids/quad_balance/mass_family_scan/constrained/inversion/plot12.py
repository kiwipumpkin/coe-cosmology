import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
R=json.load(open('step12_results.json'));from obs import ClP
l=np.arange(31);L=l[2:]
for k in['red','REF']:
    r=R[k];f,ax=plt.subplots(figsize=(9,5.5))
    ax.plot(L,np.array(r['a_motion'])[2:],'o-',label='motion / recombination path term')
    ax.plot(L,np.array(r['a_grav'])[2:],'s-',label='gravity (tidal kernel) term, A_g=1')
    ax.plot(L,np.array(r['a_sum'])[2:],'k--',label='sum')
    b=r['full']['best_axis_Ag1'];n=r['full']['nominal_axis']
    ax.plot(L,np.array(b['a_obs'])[2:],'r^-',label='Planck SMICA a_l0, best axis (l,b)=(%.1f,%.1f), chi2=%.1f'%(*b['lb'],b['chi2']))
    ax.plot(L,np.array(n['a_obs'])[2:],'mv:',label='Planck SMICA a_l0, nominal axis (0,15), chi2=%.1f'%n['chi2_Ag1'])
    ax.fill_between(L,-np.sqrt(ClP[2:]),np.sqrt(ClP[2:]),color='gray',alpha=.2,label='±1σ cosmic variance (sqrt C_l Planck)')
    ax.set_yscale('symlog',linthresh=1e-2);ax.set_xlabel('ℓ');ax.set_ylabel('a_ℓ0 [μK], CoE-pole frame');ax.set_title(f'{k}: a_ℓ0 terms vs Planck');ax.legend(fontsize=7);f.tight_layout();f.savefig(f'al0_{k}.png',dpi=120)
