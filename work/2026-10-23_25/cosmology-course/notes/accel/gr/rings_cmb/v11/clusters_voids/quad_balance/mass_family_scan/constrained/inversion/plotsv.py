import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
R=json.load(open('jac_results.json'));f,ax=plt.subplots(1,2,figsize=(11,4.5))
for U,mk in(('red','o'),('refpoly','s')):
    ax[0].semilogy(range(1,8),R[U]['singular_values_weighted'],mk+'-',label=U);ax[1].semilogy(range(1,8),R[U]['singular_values_colnormalized'],mk+'-',label=U)
ax[0].set_title('SVD of J/σ (natural units: C_k, log10 frac, DRV, A_g)');ax[1].set_title('SVD of J/σ with unit-norm columns');[a.set_xlabel('index') or a.legend() for a in ax];f.tight_layout();f.savefig('singular_values.png',dpi=120)
