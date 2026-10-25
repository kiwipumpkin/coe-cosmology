CC=[-0.041,-0.035,0.1353,0.0];Rb=1.5*(1+1/2.937)
shp=lambda x:np.exp(sum(c*(np.minimum(np.abs(x),Rb)-1)**(k+1) for k,c in enumerate(CC)))
