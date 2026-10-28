# Splits the first-draft combined image into its two original halves, unchanged.
from PIL import Image
im=Image.open('chris_fig8_orientation.png'); w,h=im.size; x=3043
im.crop((0,0,x,h)).save('fig8_orientation_seq.png'); im.crop((x,0,w,h)).save('fig8_migkas_maps.png')
