import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
R=json.load(open('signed_multipoles.json'));g=lambda k:np.array(R[k]['a'])
tot=g('b_1100_k0.npz');kin=g('none_1100_k0.npz');grav=tot-kin;bal=g([k for k in R if k.endswith('k0.45.npz')][0])
P=json.load(open('../reproduce_cmb/multipoles_vs_planck.json'));pl=np.array([r['planck_1e5'] for r in P])
l=np.arange(1,13);w=.18;fig,ax=plt.subplots(figsize=(15,7))
for i,(y,lab,c) in enumerate(((grav,'gravitational (mass − no-mass), k=0','C3'),(kin,'kinematic/emission-time only (no mass), k=0','C0'),(tot,'net, k=0 (CatWISE baseline)','k'),(bal,'net, k=0.45 (l=2 balanced)','C2'),(pl,'Planck 2018 (per-l rms)','gold'))):
    ax.bar(l+(i-2)*w,np.abs(y[1:] if len(y)==13 else y),w,label=lab,color=c)
ax.set_yscale('log');ax.set_xticks(l);ax.set_xlabel('multipole l');ax.set_ylabel('|amplitude| (1e-5, per-l rms; ×27.3 → μK)');ax.legend()
ax.set_title('CMB multipoles: gravitational vs kinematic split (dipole shown before motion removal)')
plt.savefig('quad_balance_bars.png',dpi=110,bbox_inches='tight')
print('grav',grav[1:4],'kin',kin[1:4])
