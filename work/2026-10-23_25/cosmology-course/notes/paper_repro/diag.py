exec(open('cmbann.py').read().split("meanT=lambda")[0])
for Rmpc in (1,10,100,2000):
    Rr=Rmpc/1e3/DG; R=rays(Rr,90)
    for r in R[:2]:
        qe=q_of_r(Rr,r['te']); print(Rmpc,r['beta'],'te',r['te'],'qe',qe,lay(qe,r['te']),lay(qe,1.0),'S',r['S'],'Dt(te)',Dt(r['te']))
