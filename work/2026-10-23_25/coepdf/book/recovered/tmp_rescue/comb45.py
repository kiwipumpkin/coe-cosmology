import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
top=Image.open('combined_v5.png').crop((0,70,2006,1175))
top.save('/tmp/top.png')
rows=[('A','CMB dipole cold','84, −48','cold','< 0','+0.02','✗','Our motion → distance/velocity'),
('A','Bulk flow away','107, −8','cold','< 0','−2.27','✓','Recession velocity → distance'),
('A','CMB Cold Spot','209, −57','cold','< 0','−2.33','✓','Eridanus supervoid → distance'),
('A','CMB dipole hot','264, +48','hot','> 0','+1.54','✓','Our motion → distance/velocity'),
('A','Bulk flow toward','287, +8','hot','> 0','+1.90','✓','Recession velocity → distance'),
('B','NVSS cold','16, −46','cold','< 0','−4.87','✓','Source density → sample'),
('B','Hemi. asym. cold','50, +20','cold','< 0','−2.11','✓','Power direction → sample'),
('B','CatWISE cold','58, −31','cold','< 0','−4.87','✓','Matter dipole → sample'),
('B','NVSS hot','196, +46','hot','> 0','+2.03','✓','Source density → sample'),
('B','Hemi. asym. hot','230, −20','hot','> 0','+1.62','✓','Power direction → sample'),
('B','CatWISE hot','238, +31','hot','> 0','+1.39','✓','Matter dipole → sample'),
('B','Loop I / NPS','329, +17.5','hot','> 0','+1.60','✓','Soft X-ray foreground → slope')]
fig,ax=plt.subplots(figsize=(18,5.2));ax.axis('off')
t=ax.table(cellText=[r[:7]+(r[7],) for r in rows],colLabels=['Tied to','Spot','ℓ, b','Hot/cold','Needs','Residual at spot (%)','Match','Why A or B'],loc='center',cellLoc='center',colWidths=[.06,.15,.1,.07,.06,.13,.06,.3])
t.auto_set_font_size(False);t.set_fontsize(12);t.scale(1,1.6)
for i,r in enumerate(rows,1):
  t[i,3].get_text().set_color('red' if r[3]=='hot' else 'blue')
  t[i,6].get_text().set_color('green' if r[6]=='✓' else 'red')
ax.set_title('Axis lifted 15°, then 45° counterclockwise about CoE axis  —  Score: A 4/5, B 7/7, total 11/12',fontsize=15)
fig.savefig('/tmp/tab.png',dpi=110,facecolor='white',bbox_inches='tight')
ims=[Image.open(p).convert('RGB') for p in ['/tmp/top.png','overlay_lift15_r45.png','/tmp/tab.png']]
W=max(i.width for i in ims)
ims=[i.resize((W,int(i.height*W/i.width))) for i in ims]
out=Image.new('RGB',(W,sum(i.height for i in ims)),'white');y=0
for i in ims: out.paste(i,(0,y));y+=i.height
out.save('combined_r45.png')
