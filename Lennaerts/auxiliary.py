# Auxiliary function that do not have a home elsewhere.
import numpy as np

# Set the grid in space and angle.
def setGrid(pni):
    xs = np.linspace(-np.pi,np.pi,pni[1],endpoint=False)
    phis = np.linspace(-np.pi,np.pi,pni[2],endpoint=False)
    return xs, phis

# Set the arrays that are used in the ETD1 time-stepping scheme.
def setAB(pn,phis,pni):
    n = pni[1]
    m = pni[2]
    h = pn[0]
    A = np.empty((n,n,m),dtype=np.cdouble)
    B = np.empty((n,n,m),dtype=np.cdouble)
    w = np.concatenate((np.linspace(0.0,n/2,int(n/2+1)),np.linspace(-n/2+1,-1,int(n/2-1))))
    w[int(n/2)] = 0.0 # Note that the wave number n/2 component of the first derivative must be set to zero.
    print('h is %e' % (h))
    for i1 in range(n):
        for i2 in range(n):
            Q = - (w[i1] * np.cos(phis) + w[i2] * np.sin(phis))
            A[i1,i2,:] = np.exp(h* 1J * Q)
            B[i1,i2,:] = h * np.exp(1J * Q * h/2) * np.sinc(Q*h/(2.0*np.pi))  # Using the sinc function to avoid division by zero! Note that np.sinc(x)=sin(pi x)/(pi x)                
    return A,B

# Approximate the integral int(S[T] u d phi) with the trapezoidal rule.
def quadSu(n,m,S,u):
    R = np.empty((n,n,m))
    for i1 in range(n):
        for i2 in range(n):
            for j2 in range(m):
                R[i1,i2,j2] = np.dot(S[i1,i2,:,j2],u[i1,i2,:])*2.0*np.pi/float(m)
    return R

# Approximate the integral int(S[T(x1,x2,phi,phi')] d phi'):
def quadS(n,m,S):
    Lambda = np.empty((n,n,m))
    for i1 in range(n):
        for i2 in range(n):
            for j1 in range(m):
                Lambda[i1,i2,j1] = np.sum(S[i1,i2,j1,:])*2.0*np.pi/float(m)
    return Lambda

# Approximation to the periodic Gaussian (truncated at 2*tru+1 terms):
def g(theta,sg,pni):
    M = pni[3]
    pg = 0.0
    for z in range(-M,M+1):
        pg += np.exp(-((theta+2.0*np.pi*float(z))/sg)**2)
    pg /= sg*np.sqrt(np.pi)
    return pg
