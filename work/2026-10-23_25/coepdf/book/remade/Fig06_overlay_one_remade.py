# REMADE (not original): stack of the three Migkas Fig 8 overlays (+15 tilt / +gradient / +80.8 roll) made by
# figs_letter_exact_python makers (Fig06_overlay_matched*.png), with row labels and one shared title.
from PIL import Image,ImageDraw,ImageFont
P='/workspace/figs_letter_exact_python/'
fs=['Fig06_overlay_matched.png','Fig06_overlay_matched_grad.png','Fig06_overlay_matched_grad_rot81.png'];lab=['+15° tilt','+15° tilt + gradient','+15° tilt + gradient + 80.8° roll']
ims=[Image.open(P+f).convert('RGB') for f in fs];W=max(i.width for i in ims);ims=[i.resize((W,int(i.height*W/i.width))) for i in ims]
try:F=ImageFont.truetype('DejaVuSans-Bold.ttf',26);T=ImageFont.truetype('DejaVuSans.ttf',32)
except Exception:F=T=ImageFont.load_default()
out=Image.new('RGB',(W,70+sum(i.height+40 for i in ims)),'white');d=ImageDraw.Draw(out);d.text((W//2,30),'Migkas Fig. 8 — θ=75°  A/A_all (their numbered colorbar)',font=T,fill='black',anchor='mm');y=70
for i,l in zip(ims,lab):d.text((10,y+5),l,font=F,fill='black');out.paste(i,(0,y+40));y+=i.height+40
out.save('Fig06_overlay_one_remade.png')
