import pymupdf as f,sys,re
def lines(p):
  L=[]
  for b in p.get_text('dict')['blocks']:
    for l in b.get('lines',[]):
      t=''.join(s['text'] for s in l['spans'])
      if t.strip(): L.append((round(l['bbox'][1],1),round(l['bbox'][0],1),t))
  return sorted(L)
def imgs(p): return sorted(tuple(round(x,1) for x in i['bbox']) for i in p.get_image_info())
N=lambda t:re.sub(r'\b[AWS]?\d+\b','#',t)
def compare(pa,pb):
  iss=[]
  if pa.rect!=pb.rect: iss.append('rect')
  if imgs(pa)!=imgs(pb): iss.append('images moved')
  A=lines(pa);B=lines(pb)
  if len(A)!=len(B): iss.append(f'linecount {len(A)}->{len(B)}');return iss,0
  ch=0
  for a,b in zip(A,B):
    if abs(a[0]-b[0])>0.2: iss.append('y '+a[2][:40]);break
    if N(a[2])!=N(b[2]): iss.append('text '+a[2][:50]+' || '+b[2][:50])
    if a[2]!=b[2]: ch+=1
  return iss,ch
if __name__=='__main__':
  a=f.open(sys.argv[1]);b=f.open(sys.argv[2]);offa=int(sys.argv[3]) if len(sys.argv)>3 else 0
  bad=0
  for i in range(len(b)):
    iss,ch=compare(a[i+offa],b[i])
    if iss: bad+=1;print('PAGE',i+1,iss[:4])
  print('pages',len(b),'bad',bad)
