# All parameters are set here before the time-stepping commences.
import numpy as np
from scipy.special import erf # Appears in the norms of the kernels.

def setParameters():
    # Dimensional system parameters:
    gamma = 1.0     # Speed of the organism's motion without interactions.
    L = 2.0 * np.pi # System size.
    qa = 0.0        # Relative amplitude, width and distance for the distance kernels.
    qr = 1.0
    ql = 0.0
    ma = 0.5
    mr = 0.1
    ml = 0.3
    da = 0.5
    dr = 0.0
    dl = 0.4
    sat = 0         # 0=Fetecau model, 1=Fetecau model with saturation - if sat=0 then mu, sig and mx are not referenced.
    mu = 0.3        # Threshold value (mu) and width (sig) and maximum (mx) of the sigmoidal function.
    sig = 0.1
    mx = 1.5
    ka = 1.0        # Parameter that appears in the propability functions w_a, w_r and w_l.
    kr = -1.0
    kl = 1.0
    sg = 1.0        # Width of the approximation to the delta function in the weight function.
    # Total mass - should be conserved:
    mass = 10000.0
    # Initial condition:
    ic = 0          # (currently, 0=single Gaussian bump, 1=two Guassian bumps, 2=random but smooth, change in setInitialConditions.py)
    # Visualization options:
    arrows = 1      # 0=no arrows, 1=arrow drawn at (x1,x2) for each local max phi* of u(x1,x2,phi) with length in proportion to u(x1,x2,phi*).
    # Echo parameters in dimensional units:
    print('System parameters are gamma=%f, mass=%f' % (gamma,mass))
    print('Parameters of the distance kernels:')
    print('qa=%f qr=%f ql=%f' % (qa,qr,ql))
    print('ma=%f mr=%f ml=%f' % (ma,mr,ml))
    print('da=%f dr=%f dl=%f' % (da,dr,dl))
    print('Parameters of the orientation kernel:')
    print('ka=%f kr=%f kl=%f sig=%f' % (ka,kr,kl,sg))
    if sat:
        print('Parameters of the sigmoidal function:')
        print('mu=%f sig=%f max=%f' % (mu,sig,mx))
    else:
        print('Saturation is turned off - running the Fetecau model.')
    print('mass=%e' % (mass))
    
    # Normalize using gamma (length/time) and L/2 pi (length):
    l_scale = L/(2.0*np.pi)
    t_scale = l_scale/gamma
    qa *= t_scale
    qr *= t_scale
    ql *= t_scale
    ma /= l_scale
    mr /= l_scale
    ml /= l_scale
    da /= l_scale
    dr /= l_scale
    dl /= l_scale
    mu /= l_scale
    sig /= l_scale
    mass *= l_scale**2
    
    # Auxiliary parameters (norm of distance kernel):
    nrma = np.pi * ma * (ma*np.exp(-da**2/ma**2) + np.sqrt(np.pi) * da * (1.0 + erf(da/ma)))
    nrmr = np.pi * mr * (mr*np.exp(-dr**2/mr**2) + np.sqrt(np.pi) * dr * (1.0 + erf(dr/mr)))
    nrml = np.pi * ml * (ml*np.exp(-dl**2/ml**2) + np.sqrt(np.pi) * dl * (1.0 + erf(dl/ml)))
    
    # Time-stepping parameters:
    T = 2*L/gamma # Setting T= X L/gamma means to integrate as long as it takes to cross the domain T times at speed gamma.
    Delta = 0.05  # Time step size. Instability may occur when set too large.
    nsteps = int(T/Delta) + 1 # Number of time steps.
    h = T/float(nsteps)       # Time step smaller than or equal to Delta, set so that nstep2 * h = T. 
    print('Time-stepping parameters are T=%f and h=%f (%d steps).' % (T,h,nsteps))
    # Normalize:
    h /= t_scale
    
    # Truncation parameters:
    n = 64 # spatial truncation
    m = 32 # angular truncation. These can be different but are best set to an integer power of 2 for efficient FFTs.
    dx = L/float(n)
    dphi = 2.0*np.pi/float(m)
    tru = 5 # Truncation of the series that defines the periodic Gaussian. A number greater than this seems overkill.
    print('Setting up a %d X %d point grid in space and a %d point grid in the angle...' % (n,n,m))
    print('... grid spacing is %e in space and %f in the angle.' % (dx,dphi)) 
    dx /= l_scale
    
    ps = np.array([qa,qr,ql,ma,mr,ml,da,dr,dl,ka,kr,kl,sg,mu,sig,mx,mass])
    pn = np.array([h,dx,dphi])
    pni = np.array([nsteps,n,m,tru,arrows, sat,ic])
    pa = np.array([nrma,nrmr,nrml])

    return ps, pn, pni, pa
