#input the libraries
from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from grid_setup import *
from fouriergrid_setup import *
from newtimestep import *

# Set up initial data:
N = 64
m = 64
dt = 0.05
g1=1
g2=1
G=g1+0.5*g2
gamma=0.5

#Set up spacial grid
[xx,yy]=grid_setup(N)

#initial consition with constant angle
u1=[exp(cos(pi*yy))]*m
u_hat=fft.fft2(u1) #2D fourier transform 
###Do I need to shift in this case? 
##Alt thought 
u1_hat=fft.fft(exp(cos(pi*yy)))
u_hat=[u1_hat]*m

#Probability function 
mu, sigma = 0, 0.5
w=random.normal(mu, sigma, m)
F_list=-G+G*2*pi/N*w

#set up fourier grid
[kk, ll, plus, minus]=fouriergrid_setup(N)

#used for plot
udata1 = [u1]

##Loop starts
tmax =1
nplt = floor((tmax/0.01)/dt)
nmax = int(round(tmax/dt))

for t in range(0,nmax):
    for q in range(0, m-1):
        F=F_list[q]
        p=q-1
        r=q+1
        if q==0:
            p=m-1
        if q==m-1:
            r=1
        uq_hat=u_hat[q]
        ur_hat=u_hat[r]
        up_hat=u_hat[p]
        unew_hat=newtimestep(uq_hat, ur_hat, up_hat, F, gamma, dt, plus, minus)
        u_hat[q]=unew_hat
        u1[q] = real(fft.ifft(u_hat[q]))

        


