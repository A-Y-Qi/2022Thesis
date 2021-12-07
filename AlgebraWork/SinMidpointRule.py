from numpy import *
from MidpointRule import *
import matplotlib.pyplot as plt


N = 256
k = (2*pi/N)*arange(-N/2.0,N/2.0)

fourier = lambda y: 2*cos(y)*exp(-2*pi*1j*k*y)

F_hat=midpoint(fourier, -pi, pi, 64)

plt.plot(k, real(F_hat))



