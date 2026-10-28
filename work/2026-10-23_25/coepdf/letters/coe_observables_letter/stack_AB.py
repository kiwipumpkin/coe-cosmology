# Stacks the first-draft A|B image vertically (A on top), with no other changes.
from PIL import Image; import numpy as np
im=Image.open('chris_AB_residual_spots.png'); w,h=im.size
a=np.asarray(im.convert('RGBA')); ink=((a[...,:3].min(2)<245)&(a[...,3]>0)).sum(0)
c=w//2; x=min(range(int(w*.4),int(w*.6)),key=lambda i:(ink[i],abs(i-c)))
A=im.crop((0,0,x,h)); B=im.crop((x,0,w,h))
out=Image.new('RGBA',(max(A.width,B.width),2*h),(255,255,255,0)); out.paste(A,(0,0)); out.paste(B,(0,h)); out.save('chris_AB_stacked.png'); print(x,ink[x])
