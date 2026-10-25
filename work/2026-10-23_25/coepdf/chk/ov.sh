cd /workspace/coepdf; rm -f chk/a-*; pdftoppm -r 40 -png main_rev5.pdf chk/a; python3 -c "
from PIL import Image;import glob
f=sorted(glob.glob('chk/a-*.png'));ims=[Image.open(x) for x in f];w,h=ims[0].size
for k in range(0,len(ims),5):
  S=Image.new('RGB',(w*5,h),'white')
  for i,im in enumerate(ims[k:k+5]):S.paste(im,(i*w,0))
  S.save(f'chk/ov_'+__import__('os').environ.get('T','x')+f'_{k//5}.png')"
