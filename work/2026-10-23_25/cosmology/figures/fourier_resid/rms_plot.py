exec(open('cum2.py').read().split("cols=")[0])
R={}
for nm,Y,sY,base,k in (('A',Araw,sA,AF,0),('B',Bv,sB,BG,1)):
  m=base.copy();v=[np.sqrt(np.mean((Y-1)**2))]
  for i in range(5):
    v.append(np.sqrt(np.mean((Y-m)**2)))
    if i<4: m=m+term(steps[i],th,k,1.)
  R[nm]=(v,np.mean(sY))
plt.rcParams.update({'font.size':16})
fig,ax=plt.subplots(figsize=(12,6.5),facecolor='white')
x=np.arange(6);lab=['flat (=1)','base model','+n=1','+n=2','+n=3','+n=4']
for nm,c in (('A','tab:blue'),('B','tab:red')):
  v,s=R[nm];ax.plot(x,v,'o-',color=c,lw=3,ms=10,label=f'{nm} residual RMS')
  ax.axhline(s,color=c,ls=':',lw=2,label=f'{nm} mean σ ({s:.3f})')
  for xi,vi in zip(x,v):ax.annotate(f'{vi:.3f}',(xi,vi),textcoords='offset points',xytext=(0,10 if nm=='A' else -22),ha='center',color=c,fontsize=13)
ax.set_xticks(x);ax.set_xticklabels(lab);ax.set_ylabel('RMS of (data − model)');ax.set_ylim(0,0.11)
ax.set_title('Residual RMS vs angle: A/A_all and B/B_all, 71 Migkas bins')
ax.legend(loc='upper center',bbox_to_anchor=(.5,-.1),ncol=4,frameon=False,fontsize=13)
fig.tight_layout();fig.savefig('rms_AB.png',dpi=130,facecolor='white')
