import re,sys
def arg(s,i):
    # s[i]=='{'
    d=0
    for j in range(i,len(s)):
        if s[j]=='{':d+=1
        elif s[j]=='}':
            d-=1
            if d==0:return s[i+1:j],j+1
def captions(path):
    s=open(path).read();out={}
    for m in re.finditer(r'\\(FIGc?)\{',s):
        i=m.end()-1;a,i=arg(s,i);b,i=arg(s,i);c,i=arg(s,i);out.setdefault(a,c)
    for m in re.finditer(r'\\includegraphics(\[[^\]]*\])?\{([^}]*)\}',s):
        f=m.group(2).split('/')[-1]
        k=s.find('\\caption{',m.end())
        if k>0 and k-m.end()<3000:
            c,_=arg(s,k+8);out.setdefault(f,c)
    for m in re.finditer(r'\\fig\{[^}]*\}\{[^}]*\}\{([^}]*)\}',s):
        f=m.group(1).split('/')[-1]
        k=s.find('\\caption{',m.end())
        if k>0 and k-m.end()<3000:
            c,_=arg(s,k+8);out.setdefault(f,c)
    return out
if __name__=='__main__':
    for k,v in captions(sys.argv[1]).items():print('##',k,'::',v[:150])
