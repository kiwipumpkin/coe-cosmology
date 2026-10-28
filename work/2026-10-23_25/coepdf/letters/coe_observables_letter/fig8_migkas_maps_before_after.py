# Side-by-side: current letter fig8_migkas_maps.png (left) vs fig8_migkas_maps_clean.png (right).
from PIL import Image, ImageDraw, ImageFont
A=Image.open('fig8_migkas_maps.png').convert('RGB'); B=Image.open('fig8_migkas_maps_clean.png').convert('RGB')
F=ImageFont.truetype('DejaVuSans-Bold.ttf',48)
C=Image.new('RGB',(A.width+B.width+60,max(A.height,B.height)+90),'white'); d=ImageDraw.Draw(C)
d.text((A.width//2,45),'BEFORE (letter fig. 5)',font=F,fill='black',anchor='mm'); d.text((A.width+60+B.width//2,45),'AFTER (_clean)',font=F,fill='black',anchor='mm')
C.paste(A,(0,90)); C.paste(B,(A.width+60,90)); d.line([(A.width+30,90),(A.width+30,C.height)],fill='gray',width=3)
C.save('fig8_migkas_maps_before_after.png')
