# Fourier to real and real to Fourier for space and angle. Functionality is different depending on the shape of he input matrix.
import numpy as np

small = 1e-11 # If a Fourier -> real transform gives a residual imaginary part greater than this a warning is displayed.

# Alternating 1 and -1. Elementwise product of with this array corresponds to a shift over half the domain.
def set_mask(n):
    mask = np.empty((n,n))
    for i1 in range(n):
        for i2 in range(n):
            mask[i1,i2] = (-1)**(i1+i2)
    return mask

# Normalized and shifted so that fh(k) is the coefficient of orthogonal projection of f(x) onto exp(i k x) on [-pi,pi).
def r2f_space(n,m,u):
    mask = set_mask(n)
    uh = np.empty_like(u,dtype=np.cdouble)
    if np.ndim(u) == 2:
        uh[:,:] = np.fft.fft2(u[:,:]) / float(n**2)   # First DFT...
        uh[:,:] *= mask                               # then shift because the grid spans [-pi,pi)
    elif np.ndim(u) == 3:
        for j in range(m):
            uh[:,:,j] = np.fft.fft2(u[:,:,j]) / float(n**2)
            uh[:,:,j] *= mask
    elif np.ndim(u) == 4:
        for j1 in range(m):
            for j2 in range(m):
               uh[:,:,j1,j2] = np.fft.fft2(u[:,:,j1,j2]) / float(n**2)
               uh[:,:,j1,j2] *= mask
    return uh

def f2r_space(n,m,uh):
    mask = set_mask(n)
    u = np.empty_like(uh,dtype=np.double)
    dum = np.empty((n,n),dtype=np.cdouble)
    if np.ndim(uh) == 2:
        dum = np.fft.ifft2(uh[:,:] * mask) * float(n**2)  # Inverse of r2f: first shift then IFFT
        cmplx = np.linalg.norm(np.imag(dum))
        if cmplx > small:
            print('Imaginary part of %e in f2r [n,n].' % (cmplx))
        u = np.real(dum)
    elif np.ndim(uh) == 3:
        cmplx = 0.0
        for j in range(m):
            dum = np.fft.ifft2(uh[:,:,j] * mask) * float(n**2) 
            cmplx = np.max([cmplx,np.linalg.norm(np.imag(dum))])
            u[:,:,j] = np.real(dum)
        if cmplx > small:
            print('Imaginary part of %e in f2r [n,n,m].' % (cmplx))

    elif np.ndim(uh) == 4:
        cmplx = 0.0
        for j1 in range(m):
            for j2 in range(m):
                dum = np.fft.ifft2(uh[:,:,j1,j2] * mask) * float(n**2)
                cmplx = np.max([cmplx,np.linalg.norm(np.imag(dum))])
                ind = np.unravel_index(np.argmax(np.imag(dum), axis=None), dum.shape)
                u[:,:,j1,j2] = np.real(dum)
        if cmplx > small:
            print('Imaginary part of %e in f2r [n,n,m,m]; max of %e at [%d,%d,%d,%d].' % (cmplx,dum[ind[0],ind[1]],ind[0],ind[1],j1,j2))

    return u

def r2f_angle(n,m,uh):
    mask = set_mask(m)
    uhh = np.empty_like(uh,dtype=np.cdouble)
    if np.ndim(uh) == 2:
        uhh[:,:] = np.fft.fft2(uh) / float(m**2)       # First DFT...
        uhh[:,:] *= mask                               # then shift because the grid spans [-pi,pi)
    elif np.ndim(uh) == 3:
        for i1 in range(n):
            for i2 in range(n):
                uhh[i1,i2,:] = np.fft.fft(uh[i1,i2,:]) / float(m)
                uhh[i1,i2,:] *= mask[:,0]
    elif np.ndim(uh) == 4:
        for i1 in range(n):
            for i2 in range(n):
                uhh[i1,i2,:,:] = np.fft.fft2(uh[i1,i2,:,:]) / float(m**2)
                uhh[i1,i2,:,:] *= mask
    return uhh

def f2r_angle(n,m,up):
    mask = set_mask(m)
    u = np.empty_like(up,dtype=np.cdouble)
    if np.ndim(up) == 3:
        for i1 in range(n):
            for i2 in range(n):
                u[i1,i2,:]  = np.fft.ifft(up[i1,i2,:] * mask[:,0]) * float(m)
    elif np.ndim(up) == 4:
        for i1 in range(n):
            for i2 in range(n):
                u[i1,i2,:,:] = np.fft.ifft2(up[i1,i2,:,:] * mask) * float(m**2)
    return u

    
