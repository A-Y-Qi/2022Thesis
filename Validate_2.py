from numpy import *
from MidpointRule import *
import matplotlib.pyplot as plt
from scipy.linalg import dft
import scipy.special

N = 256
power=10

k = arange(0,N) #wave number 

Grid_x=zeros(N, dtype=complex)
function=zeros(N, dtype=complex)
##Grid point in physical space
for r in range (0,N):
    function[r]=1/sqrt(N)*(2*cos(2*pi*(1/(2*N)+r/N)))**power
    Grid_x[r]=2*pi*(1/(2*N)+r/N)
    
FFT=1/(sqrt(N))*fft.fft(function)
    



##From paper:
compare=arange(0,2*power+2,2)
result=zeros(N, dtype=complex)
for s in range(0,N):
    add=0
    for p in range (-N,N):
        nv_2=power-s+p*N
        #nv_2=power-s
        if nv_2 in compare:
            nv=nv_2/2
            add=add+(math.factorial(power)/(math.factorial(nv)*math.factorial(power-nv)))*exp(-2*pi*1j*1/(2*N)*(nv_2-power)) 
    result[s]=add
    

    
fig1, ax1 = plt.subplots()
ax1.plot(k, real(FFT), 'o', label='FFT')
ax1.plot(k, real(result), '.', label='Paper')
#ax.axis('equal')
leg = ax1.legend();
fig1.suptitle ("The comparison between Paper and FFT Real Part with power="+str(power))

fig2, ax2 = plt.subplots()
ax2.plot(k, imag(FFT), 'o', label='FFT')
ax2.plot(k, imag(result), '.', label='Paper')
#ax.axis('equal')
leg = ax2.legend();
fig2.suptitle ("The comparison between Paper and FFT Imag Part with power="+str(power))

fig3, ax3 = plt.subplots()
ax3.plot(k, real(FFT)-real(result), 'o', label='Real Difference')
leg = ax3.legend();
fig3.suptitle ("The difference between Paper and FFT Real Part with power="+str(power))

fig4, ax4 = plt.subplots()
ax4.plot(k, imag(FFT)-imag(result), 'o', label='Imag Difference')
leg = ax4.legend();
fig4.suptitle ("The difference between Paper and FFT Imag Part with power="+str(power))




#fig3, ax3 = plt.subplots()
#ax3.plot(Grid_x, function, 'o', label='FFT')
