import pymupdf as f,sys,numpy as np
a=f.open(sys.argv[1]);b=f.open(sys.argv[2]);print(len(a),len(b))
for i in range(min(len(a),len(b))):
  pa=a[i].get_pixmap(dpi=60);pb=b[i].get_pixmap(dpi=60)
  if (pa.w,pa.h)!=(pb.w,pb.h): print(i+1,'size');continue
  x=np.frombuffer(pa.samples,np.uint8).reshape(pa.h,pa.w,-1).astype(int);y=np.frombuffer(pb.samples,np.uint8).reshape(pb.h,pb.w,-1).astype(int)
  d=np.abs(x-y).max(2)>40;ys,xs=np.nonzero(d)
  print(i+1,d.sum(), (ys.min(),ys.max(),xs.min(),xs.max()) if d.sum() else '')
