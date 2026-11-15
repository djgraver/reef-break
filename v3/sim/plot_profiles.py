import numpy as np, sys, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
d=np.load(sys.argv[1]); t=d['t']; n=d['n']; off=np.concatenate([[0],np.cumsum(n)])
k=int(sys.argv[2]); st=int(sys.argv[3]); skip=int(sys.argv[4]) if len(sys.argv)>4 else 0
last=len(t)-1-skip
fig,ax=plt.subplots(1,2,figsize=(14,6))
for i in range(max(0,last-k*st),last+1,st):
    x=d['x'][off[i]:off[i+1]]; y=d['y'][off[i]:off[i+1]]
    ax[0].plot(x,y,'-',lw=.8); 
x=d['x'][off[last]:off[last+1]]; y=d['y'][off[last]:off[last+1]]
ax[1].plot(x,y,'-',lw=.6,marker='.',ms=2.5)
j=np.argmax(y); xc=x[j]
ax[0].plot(d['bx'],d['by'],'k'); ax[0].set_xlim(xc-14,xc+12); ax[0].set_ylim(-2.5,4); ax[0].set_aspect('equal'); ax[0].grid(alpha=.3)
ax[1].set_xlim(xc+1,xc+8); ax[1].set_ylim(y.max()-4.2,y.max()+.3); ax[1].set_aspect('equal'); ax[1].grid(alpha=.3); ax[1].set_title('t=%.2f N=%d'%(t[last],len(x)))
plt.tight_layout(); plt.savefig(sys.argv[1].replace('.npz','_z.png'),dpi=75)
