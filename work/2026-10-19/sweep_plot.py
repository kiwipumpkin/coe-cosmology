import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
J=json.load(open('sweep/sweep_dtt.json'));PR=J['PR'];R=J['rows'];U=[float(k) for k in R];K=list(R)
fig,ax=plt.subplots(figsize=(13,8))
for t,q,col,ls,name in (('nomass','rng','C0','-','dT/T range ÷ Planck, no mass'),('mass','rng','C3','-','dT/T range ÷ Planck, with Migkas best-fit mass'),
                        ('nomass','rrng','C0','--','range left after dipole removal ÷ Planck, no mass'),('mass','rrng','C3','--','range left after dipole removal ÷ Planck, with mass')):
    y=[R[k][t][q]/PR for k in K]
    ax.plot(U,y,ls,color=col,marker='o',lw=2.5,label=f"{name}: {y[0]:.3g} at 4 c, {y[1]:.3g} at 5 c, {y[3]:.3g} at 8 c, {y[-1]:.3g} at 30 c")
ax.axhline(1,color='k',lw=.8,ls=':');ax.text(30,1.08,'Planck ring range',ha='right',fontsize=10)
ax.set_yscale('log');ax.set_xlim(0,31);ax.set_xticks(U);ax.grid(alpha=.3,which='both')
ax.set_xlabel('our radial speed (multiples of c)');ax.set_ylabel('ratio to Planck ring range (log scale)')
ax.legend(loc='lower left',bbox_to_anchor=(0,1.01),fontsize=10.5,frameon=False)
ax.set_title('green today = 4.2 Gpc at every speed (D scaled); recombination radius 0.5 D; mass = Migkas best fit at each speed (0.0166–0.0170 of critical)',fontsize=9.5,pad=-14)
fig.savefig('dTT_vs_speed_c4sweep.png',dpi=110,bbox_inches='tight')
