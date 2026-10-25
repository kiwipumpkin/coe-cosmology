tectonic main_rev5.tex 2>&1 | grep -i "error" ; pdfinfo main_rev5.pdf|grep Pages; rm -rf rev5_pages; mkdir rev5_pages; pdftoppm -r 45 -png main_rev5.pdf rev5_pages/p; python3 -c "
from PIL import Image;import glob
fs=sorted(glob.glob('rev5_pages/p-*.png'));ims=[Image.open(f) for f in fs];w,h=ims[0].size
n=len(ims);c=5;r=(n+c-1)//c;S=Image.new('RGB',(w*c,h*r),'white')
for i,im in enumerate(ims):S.paste(im,((i%c)*w,(i//c)*h))
S.save('chk/sheet_$1.png')"
