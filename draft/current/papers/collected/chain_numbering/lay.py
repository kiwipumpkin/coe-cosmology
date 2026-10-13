import pymupdf as f,sys,re
def words(p): return [(round(w[0],1),round(w[1],1),round(w[2],1),round(w[3],1),w[4]) for w in p.get_text('words')]
def imgs(p): return sorted(tuple(round(x,1) for x in i['bbox']) for i in p.get_image_info())
def compare(pa,pb,tag):
  issues=[];chg=[]
  if pa.rect!=pb.rect: issues.append('rect')
  if imgs(pa)!=imgs(pb): issues.append('images moved')
  A=words(pa);B=words(pb)
  sa=set(A);sb=set(B)
  da=[w for w in A if w not in sb]; db=[w for w in B if w not in sa]
  # lines: compare y coordinates of lines
  la=sorted(set(w[1] for w in A)); lb=sorted(set(w[1] for w in B))
  if la!=lb: issues.append('line y changed')
  for w in da+db:
    if not re.search(r'\d',w[4]): issues.append('nonnum word changed:'+w[4])
  nonnum=[w for w in da if not re.search(r'\d',w[4])]
  return issues,[w[4] for w in da],[w[4] for w in db]
if __name__=='__main__':
  a=f.open(sys.argv[1]);b=f.open(sys.argv[2]);offa=int(sys.argv[3]) if len(sys.argv)>3 else 0
  for i in range(len(b)):
    iss,x,y=compare(a[i+offa],b[i],i)
    if iss or x: print(i+1,iss[:6],'OLD',x[:12],'NEW',y[:12])
