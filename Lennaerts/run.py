# Code for the simulation of Fetecau's 2+1D aggregation model with saturation.
# Equations scales with length scale L/2 pi and time scale L/(2 pi gamma).
# The emphasis is on transparency and correctness, not on efficiency.
import numpy as np
from setParameters import setParameters               # This is where all parameters are set.
from setKernels import setKernels                     # This is where the kernels are pre-computed.
from auxiliary import *                               # Some auxiliary functions like approximations for integrals.
from setInitialCondition import setInitialCondition   # Where the initial condition is set.
from march import march                               # The actual time-stepping.
from plotSolution import plotSolution                 # Makes a density plot with arrows based on u(x1,x2,phi).
import matplotlib.pyplot as plt

# First set all the parameters, returned in Parameters of the System (ps), Parameters of the Numerics (pn),
# Parameters of the Numerics that are Integers (pni) and Parameters that are Auxiliary (like the normalization of the periodic Gaussian).
ps, pn, pni, pa = setParameters()
# Order of ps is qa, qr, ql, ma, mr, ml, da, dr, dl, ka, kr, kl, sig, mx, mass.
# Order of pn is h, dx, dphi.
# Order of pni is nsteps, n, m, tru, arrows, sat, ic.
# Order of pa is nrma, nrmr, nrml.
# See setParameters.py for a description of each parameter.

# Grid in space and angle:
xs, phis = setGrid(pni)

# Set auxiliary variables for the ETD1 method:
A,B = setAB(pn,phis,pni)

# Compute the kernels:
Qa, Qr, Ql = setKernels(ps,pni,pa,xs,phis)

# Initialize the solution:
u = setInitialCondition(ps,pn,pni,xs,phis,pni[6])

# Now time-step:
out = march(ps,pn,pni,u,A,B,Qa,Qr,Ql)

# Write output to file for off-line plotting:
File = open("timeseries","w")
np.savetxt(File,out,fmt='%16.9e',delimiter=' ')

