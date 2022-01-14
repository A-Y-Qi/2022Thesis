from numpy import *
from MidpointRule import *
import matplotlib.pyplot as plt
from scipy.linalg import dft
import scipy.special


N = 7000
N2=7000


k = arange(-N/2.0,N/2.0) #wave number 
k2=arange(-N2/2.0,N2/2.0) #wave number for CFT
m = (2*pi/N)*arange(-N/2.0,N/2.0) #grid number
function=exp(cos(m))
h=2*pi/N
##function itself y=exp(cos(x))
##power=1

##analytic solution 
fourier = lambda x: exp(cos(x))*exp(-1*1j*k2*x)
##this should be the anayticak answer
F_hat=midpoint(fourier, -pi, pi, N2 ) #CFT result

##Find the DFT answer 
DFT=zeros(N, dtype=complex)
for i in range(0,N-1):
    w_number = k[i]
    sum_n=0
    for j in range(0,N-1):
        sum_n=sum_n+exp(-1j*w_number*m[j])*function[j]
    DFT[i]=h*sum_n
    
start=k[0]
end=k[-1]
index_0=where(k2==start)
index_end=where(k2==end)
index_0=index_0[0][0]
index_end=index_end[0][0]
diffrc_at_end=sum(F_hat[0:index_0])+sum(F_hat[index_end:-1])
temp=F_hat[index_0:index_end+1]
diffrc_at_point_1=real(F_hat[index_0:index_end+1])-real(DFT)
diffrc_at_point_2=imag(F_hat[index_0:index_end+1])-imag(DFT)


        




#plt.plot(k2, real(F_hat))
#plt.plot(k, real(DFT),".")
#plt.plot(k, real(F_hat-DFT),"r.")

fig1, ax1 = plt.subplots()
ax1.plot(k2, real(F_hat), 'o', label='cutted CFT with N=7000')
ax1.plot(k, real(DFT), '.', label='DFT with N =1024')
#ax.axis('equal')
leg = ax1.legend();
plt.xlim(-N, N)
fig1.suptitle ("The comparison between CFT and DFT Real Part")

fig2, ax2 = plt.subplots()
ax2.plot(k2, imag(F_hat), 'o', label='cutted CFT with N=7000')
ax2.plot(k, imag(DFT), '.', label='DFT with N =1024')
leg = ax2.legend();
plt.xlim(-N, N)
fig2.suptitle ("The comparison between CFT and DFT Img Part")

fig3, ax3 = plt.subplots()
ax3.plot(k, diffrc_at_point_1, '.', label='difference at point real part')
leg = ax3.legend();
#plt.xlim(-N, N)
fig3.suptitle ("The difference at the point real part")

fig4, ax4 = plt.subplots()
ax4.plot(k, diffrc_at_point_2, '.', label='difference at point imag part')
leg = ax4.legend();
#plt.xlim(-N, N)
fig4.suptitle ("The difference at the point imag part")


