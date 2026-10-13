import json,re,sys,urllib.request,concurrent.futures as cf
pat=json.load(open('pat_refs.json'))
c={}
for s in pat:
  n=int(re.match(r'\[(\d+)\]',s).group(1))
  m=re.search(r'https://doi.org/(\S+?)(,|$)',s)
  if m: c['P%d'%n]=m.group(1)
c.update({'C1':'10.1086/300499','C2':'10.1086/307221','C3':'10.1088/0264-9381/4/6/025','C4':'10.1073/pnas.15.3.168','C5':'10.1051/0004-6361/201936602','C7':'10.1098/rspa.1950.0049','C9':'10.1073/pnas.20.3.169','C10':'10.1093/mnras/107.5-6.410','C11':'10.1051/0004-6361/201935201','C12':'10.3847/2041-8213/abdd40',
'B22':'10.1086/148982','B23':'10.1038/217511a0','P9':'10.1080/00018739400101505','P29':'10.1103/PhysRevE.76.066301','P55':'10.1088/0004-637X/759/2/78',
'N1':'10.1007/BF00759270','N2':'10.1016/S1387-6473(02)00234-8','C8r':'10.1023/A:1018855621348'})
def g(k):
  try:
    m=json.load(urllib.request.urlopen('https://api.crossref.org/works/'+urllib.request.quote(c[k]),timeout=30))['message']
    return k,c[k],(m.get('title') or [''])[0][:70],m.get('volume'),m.get('page'),m.get('issued',{}).get('date-parts',[[0]])[0][0],(m.get('author') or [{}])[0].get('family')
  except Exception as e: return k,c[k],'ERR '+str(e)
with cf.ThreadPoolExecutor(8) as ex:
  for r in ex.map(g,c): print(r)
