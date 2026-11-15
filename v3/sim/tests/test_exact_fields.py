import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np, warnings
warnings.filterwarnings('ignore')
from solver import *
L=100.; N=128; x=np.arange(N)*L/N; k=2*np.pi/L
reef=Reef(L, lambda x: 10-7*np.exp(-((x-60)/8)**2))
sol=Solver(reef); h=reef.h
for nm in (1,3,8,20,35):
    for zeta in (x+0j, x+1j*0.5*np.cos(k*x), reef.inv(x+0j)):
        W=np.cos(nm*k*(zeta+1j*h))/np.cosh(nm*k*h); Wz_ex=-nm*k*np.sin(nm*k*(zeta+1j*h))/np.cosh(nm*k*h)
        dz,dp,Wz,ze,fp,z=sol.derivs(zeta,W.real)
        print(nm, "rel err %.2e"%(np.abs(Wz-Wz_ex).max()/np.abs(Wz_ex).max()), end=' | ')
    print()
