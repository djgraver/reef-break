import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np, warnings
from solver import *
L=100.; h0=10.
N=128; x=np.arange(N)*L/N; k=2*np.pi/L
a=0.05; w=np.sqrt(g*k*np.tanh(k*h0)); T=2*np.pi/w; dt=T/200
eta=a*np.cos(k*x); phi0=a*g/w*np.sin(k*x)
def go(reef,nst,a_scale=1):
    sol=Solver(reef); zeta=reef.inv(x+1j*eta*a_scale); phi=phi0*a_scale
    E0=sol.energy(zeta,phi)
    for i in range(nst):
        zeta,phi=sol.step(zeta,phi,dt)
        zeta=lowpass(zeta-x,36)+x; phi=lowpass(phi,36).real
    z=reef.f(zeta)
    return z,E0,sol.energy(zeta,phi)
z,E0,E1=go(Reef(L, lambda x: h0+0*x),400)
# crest location by fit of fundamental phase
ph=np.angle(np.sum(z.imag*np.exp(-1j*k*z.real)))
print("flat: 2 periods, phase err(rad)",ph,"amp",z.imag.max(),E0,E1)
z,E0,E1=go(Reef(L, lambda x: 10-7*np.exp(-((x-60)/8)**2)),600,4)
print("reef: ",E0,E1, z.imag.max(), z.imag.min())
