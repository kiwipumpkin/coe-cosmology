exec(open("paper_repro/stretch/gr_exact.py").read().split("b_=np.loadtxt")[0])
for nm,R in out.items():
    c=np.fft.rfft(R)/len(R); a=2*np.abs(c[1:9]); ph=np.degrees(-np.angle(c[1:9]))%360
    print(nm,"RMS %.3e"%rms(R)); print("  m:amp(x1e-3)@phase "+"  ".join(f"{m}:{a[m-1]*1e3:.3f}@{ph[m-1]/m:.0f}" for m in range(1,9)))
