# Set the initial condition u(x1,x2,phi,0).
import numpy as np
from numpy.random import default_rng  # Currently the best way to use NumPy's RNGs, not backward compatible though.
from auxiliary import g               # Periodic Gaussian.
from r2f_f2r import *                 # Fourier->real->Fourier routines.

def setInitialCondition(ps,pn,pni,xs,phis,switch):
    n = pni[1]
    m = pni[2]
    mass = ps[16]
    u = np.empty((n,n,m))
    if switch == 0: # Single Gaussian bump.
        x0 = 0.4
        y0 = -0.2
        sigx = 0.7
        sigy = 0.9
        phi0 = np.pi*0.9
        sigphi = 0.5
        for i1 in range(n):
            for i2 in range(n):
                for i3 in range(m):
                    u[i1,i2,i3] = mass*g(xs[i1]-x0,sigx,pni)*g(xs[i2]-y0,sigy,pni)*g(phis[i3]-phi0,sigphi,pni)
    elif switch == 1: # Two Gaussian bumps.
        x00 = 0.4
        y00 = -0.2
        sigx0 = 0.7
        sigy0 = 0.9
        phi00 = np.pi/2.0
        sigphi0 = 0.3
        x01 = -0.4
        y01 = 0.3
        sigx1 = 0.5
        sigy1 = 0.6
        phi01 = -np.pi
        sigphi1 = 0.4
        for i1 in range(n):
            for i2 in range(n):
                for i3 in range(m):
                    u[i1,i2,i3] = mass*g(xs[i1]-x00,sigx0,pni)*g(xs[i2]-y00,sigy0,pni)*g(phis[i3]-phi00,sigphi0,pni)/2.0
                    u[i1,i2,i3] += mass*g(xs[i1]-x01,sigx1,pni)*g(xs[i2]-y01,sigy1,pni)*g(phis[i3]-phi01,sigphi1,pni)/2.0

    elif switch == 2: # Random but smooth.
        w_s = np.concatenate((np.linspace(0.0,n/2,int(n/2+1)),np.linspace(-n/2+1,-1,int(n/2-1))))
        w_a = np.concatenate((np.linspace(0.0,m/2,int(m/2+1)),np.linspace(-m/2+1,-1,int(m/2-1))))
        rng = default_rng()
        # Random, non-negative, grid point values: 
        u = rng.random((n,n,m))
        # Fourier transform in space and angle and apply an exponential mask to damp higher wave numbers:
        uh = r2f_space(n,m,u)
        uhp = r2f_angle(n,m,uh)
        # The smaller sc_s and sc_a, the smoother u(x1,x2,phi,0).
        sc_s = float(2)
        sc_a = float(2)
        for i1 in range(n):
            for i2 in range(n):
                for j1 in range(m):
                    k = np.sqrt(w_s[i1]**2+w_s[i2]**2)
                    q = np.abs(w_a[j1])
                    uhp[i1,i2,j1] *= np.exp(-k/sc_s - q/sc_a)
        sc = ps[16] / (8.0*np.pi**3 * np.real(uhp[0,0,0]))  # Mass input parameter over the integral of the random IC.
        uhp *= sc                                           # Normalize so that the total mass is as requested.
        uh = f2r_angle(n,m,uhp)
        u = f2r_space(n,m,uh)                               # Return u on the grid.
        return u
        
    else:
        u = 0.
    print('Maximal value of initial u is %f.' % (np.amax(u)))
    return u
    
