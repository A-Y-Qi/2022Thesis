# Compute the interaction terms.
import numpy as np
from r2f_f2r import *            # Fourier->real->Fourier for space and angle.
from auxiliary import *          # Varuous functions like trapezoidal integral approximations.

def S(x,mu,sig,H):
    sigm = H*(np.pi/2.0 + np.arctan((x-mu)/sig))/np.pi
    return sigm

def computeInteractions(n,m,uh,Qa,Qr,Ql,ps,pni):
    N = np.empty((n,n,m,m),dtype=np.cdouble)
    Nr = np.empty((n,n,m,m))
    # Integrate uh over the angle for copmuting attraction and repulsion:
    uph = r2f_angle(n,m,uh)
    dum = 2.0*np.pi*uph[:,:,0]
    
    # Spatial convolution of Q (product of Kd, Ko and w) and u for attraction and repulsion:
    for j1 in range(m):
        for j2 in range(m):
            N[:,:,j1,j2] = (ps[0] * Qa[:,:,j1,j2] + ps[1] * Qr[:,:,j1,j2]) * dum
    # Add alignment (this takes a few lines of algebra to derive):
    dum = np.empty_like(N)
    for i1 in range(n):
        for i2 in range(n):
            for j1 in range(m):
                for j2 in range(m):
                    s1 = np.mod(j1+j2,m)
                    s2 = np.mod(m - j2,m)
                    dum[i1,i2,j1,j2] = ps[2] * Ql[i1,i2,s1,s2] * uph[i1,i2,s1]
    N += f2r_angle(n,m,dum)
                    
    # Inverse transform to find the grid point values:
    Nr = f2r_space(n,m,N)
    Nmin = np.min(Nr)
    Nmax = np.max(Nr)

    # Evaluate S on the grid if sat=1:
    if pni[5]:
        for i1 in range(n):
            for i2 in range(n):
                for j1 in range(m):
                    for j2 in range(m):
                        Nr[i1,i2,j1,j2] = S(Nr[i1,i2,j1,j2],ps[13],ps[14],ps[15])
    
    # From S(T) now first compute the reorientation term:
    u = f2r_space(n,m,uh)
    R = quadSu(n,m,Nr,u)
    # Then integrate over the "to" angle to find the turning rate Lambda:
    Lambda = quadS(n,m,Nr)
    R = -Lambda * u + R
    Rh = r2f_space(n,m,R)

    return Rh, Nmin, Nmax
