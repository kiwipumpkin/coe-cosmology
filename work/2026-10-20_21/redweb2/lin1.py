D=np.array(d['R0.8_0.01']['Dl']);plt.figure(figsize=(10,5));plt.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='ko',ms=4,label='Planck')
plt.plot(l[2:],D[2:],'-o',color='teal',ms=3,label=f"red + web seed 2000, 0.8 Gpc, 1% (χ²={d['R0.8_0.01']['chi2']:.1f}, R²={r2(D[L]):.2f})")
plt.xlabel('ℓ');plt.ylabel('D_ℓ (μK²)');plt.xlim(1.5,30.5);plt.ylim(0);plt.legend();plt.title('Linear: Planck vs red + web seed 2000 (0.8 Gpc, 1%)');plt.tight_layout();plt.savefig('red_web_s2000_08_1pct_linear.png',dpi=120)
