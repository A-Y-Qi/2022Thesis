# Functions for pre-computing the attraction, repulsion and alignment kernels for repeated use during the time-stepping.
import numpy as np
from r2f_f2r import *    # This makes available the Fourier->real->Fourier routines.
from auxiliary import g  # Periodic Gaussian the appears in the probability function.

# Turning function:
def v(theta,ps,force):
    if force == 0:  # Attraction.
        k = ps[9]
    elif force == 1: # Repulsion.
        k = ps[10]
    elif force == 2: # Alignment.
        k = ps[11]
    return k*np.sin(theta) # Note that other forms of v are proposed in Fetecau's paper.

# Weight function:
def w(alpha,beta,ps,pni,force):
    return g(alpha - v(beta,ps,force),ps[12],pni)

# Distance kernel:        
def Kd(x1,x2,ps,pa,force):
    r = np.sqrt(x1**2+x2**2)
    if force == 0:    # attraction
        d = ps[6]
        m = ps[3]
        nrm = pa[0]
    elif force == 1:  # repulsion
        d = ps[7]
        m = ps[4]
        nrm = pa[1]
    elif force == 2:  # alignment
        d = ps[8]
        m = ps[5]
        nrm = pa[2]
    return np.exp(-(r - d)**2/m**2)/nrm

# Returns the angle between the positive x1-axis and the point (z1,z2) in the domain [-pi,pi].
def computePsi(z1,z2):
    return np.arctan2(z2,z1)

# Orientation kernel:
def Ko(z1=0,z2=0,theta=0,phi=0,ps=0,force=0):
    if force == 0:    # attraction
        psi = computePsi(z1,z2)
        K = (1.0-np.cos(phi-psi))/(2.0*np.pi)
    elif force ==1:   # repulsion
        psi = computePsi(z1,z2)
        K = (1.0+np.cos(phi-psi))/(2.0*np.pi)
    elif force ==2:   # alignment
        K = (1.0-np.cos(phi-theta))/(2.0*np.pi)
    return K

# Set the product of Kd, Ko and w:
def setKernels(ps,pni,pa,xs,phis):
    n = pni[1]
    m = pni[2]
    # First attraction:
    Qa = np.empty((n,n,m,m),dtype=np.cdouble)
    # Compute its grid point values:
    for i1 in range(n):
        for i2 in range(n):
            psi = computePsi(xs[i1],xs[i2])
            K = Kd(xs[i1],xs[i2],ps,pa,0)
            for j1 in range(m):
                for j2 in range(m):                    
                    Qa[i1,i2,j1,j2] =  K * Ko(z1=xs[i1],z2=xs[i2],phi=phis[j1],ps=ps,force=0) * w(phis[j1]-phis[j2],phis[j1]-psi,ps,pni,0)
    # Fourier transform in space only:
    Qa = r2f_space(n,m,Qa)
    # Then repulsion:
    Qr = np.empty((n,n,m,m),dtype=np.cdouble)
    for i1 in range(n):
        for i2 in range(n):
            psi = computePsi(xs[i1],xs[i2])
            K = Kd(xs[i1],xs[i2],ps,pa,1)
            for j1 in range(m):
                for j2 in range(m):                    
                    Qr[i1,i2,j1,j2] =  K * Ko(z1=xs[i1],z2=xs[i2],phi=phis[j1],ps=ps,force=1) * w(phis[j1]-phis[j2],phis[j1]-psi,ps,pni,1)
    Qr = r2f_space(n,m,Qr)
    # Then alignment. Qal is Fourier transformed in space and time since the convolution in the turning function involved both:
    Ql = np.empty((n,n,m,m),dtype=np.cdouble)
    K1 = np.empty((n,n),dtype=np.cdouble)
    for i1 in range(n):
        for i2 in range(n):
            K1[i1,i2] = Kd(xs[i1],xs[i2],ps,pa,2)    
    K1 = r2f_space(n,m,K1)
    K2 = np.empty((m,m),dtype=np.cdouble)
    for j1 in range(m):
        for j2 in range(m):
            K2[j1,j2] = Ko(phi=phis[j1],ps=ps,force=2) * w(phis[j2],phis[j1],ps,pni,2)
    K2 = r2f_angle(n,m,K2)
    # Ql is a direct product of K1 and K2. It is wasteful to construct an n X n X m X m array but that does make the code easier to parse:
    for i1 in range(n):
        for i2 in range(n):
            for j1 in range(m):
                for j2 in range(m):
                    Ql[i1,i2,j1,j2] = K1[i1,i2]*K2[j1,j2]
    return Qa, Qr, Ql
