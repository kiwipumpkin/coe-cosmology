import pymupdf as f,sys,re,lay2
N=lambda t:re.sub(r'[AWS]?\d+','#',t)
def seq(p): return N(' '.join(w[4] for w in p.get_text('words')))
def draw(p): return len(p.get_drawings())
def check(a,b,offa=0,pages=None):
  res=[]
  for i in (pages or range(len(b))):
    pa=a[i+offa];pb=b[i];r=[]
    if pa.rect!=pb.rect: r.append('RECT')
    if lay2.imgs(pa)!=lay2.imgs(pb): r.append('IMAGES MOVED')
    if sorted(seq(pa).split())!=sorted(seq(pb).split()): r.append('CONTENT DIFFERS')
    A=lay2.lines(pa);B=lay2.lines(pb)
    reflow=sum(1 for x,y in zip(A,B) if N(x[2])!=N(y[2]))+abs(len(A)-len(B))
    if reflow: r.append(f'reflow lines~{reflow}')
    if r: res.append((i,r))
  return res
if __name__=='__main__':
  a=f.open(sys.argv[1]);b=f.open(sys.argv[2])
  for i,r in check(a,b): print('page',i+1,r)
