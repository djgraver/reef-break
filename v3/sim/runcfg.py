import numpy as np, sys, warnings, json
warnings.filterwarnings('ignore')
from solver import *
S=lambda s: 0.5*(1+np.tanh(s))
def ramp(s, r=0.15):
    # linear ramp 0..1 with rounded corners (radius r)
    sm=lambda v: 0.5*(v+np.sqrt(v*v+r*r))
    return sm(s)-sm(s-1)
def make(name, L=440., h0=12., hr=1.5, x1=180., x2=285., xe=405., we=10., a=1.4, lam=140., sig=60., x0=120., rr=0.15, **kw):
    depth=lambda x: h0-(h0-hr)*ramp((x-x1)/(x2-x1),rr)*(1-S((x-xe)/we))
    def init(x):
        eta=a*np.exp(-((x-x0)/sig)**2)*np.cos(2*np.pi*(x-x0)/lam)
        d=depth(x)
        u=2*(np.sqrt(g*(d+eta))-np.sqrt(g*d))
        return eta,u
    cfg=dict(L=L, depth=depth, init=init, T=30., M=64, smooth=2.5)
    cfg.update(kw)
    return cfg
if __name__=="__main__":
    name=sys.argv[1]; P=json.loads(sys.argv[2]) if len(sys.argv)>2 else {}
    run(make(name,**P), name+'.npz')
