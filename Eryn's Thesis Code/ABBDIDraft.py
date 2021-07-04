#input the libraries
from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA

# Set up grid and two-solution initial data:
N = 64
n=64
m=64
dt = 0.05
x = (2/N)*arange(-N/2.0,N/2.0)
y = (2/N)*arange(-N/2.0,N/2.0)
h=2/N
xx,yy=meshgrid(x,y)
xx=xx.flatten()
yy=yy.flatten()
u1=exp(cos(pi*yy))
u_hat=fft.fft(u1)
uold=exp(cos(pi*yy-dt*0.2))
uold_hat=fft.fft(uold)
w=1
g1=1
g2=1
G=g1+0.5*g2
F=-G+G*2*pi/n*w
gamma=0.5
k=hstack([arange(0,N/2.0),0,arange(-N/2.0+1,0)])
l=hstack([arange(0,N/2.0),0,arange(-N/2.0+1,0)])
kk,ll=meshgrid(k,l)
kk=kk.flatten()
ll=ll.flatten()
onsarray=ones(len(kk))
udata1 = [u1]

##Loop starts
tmax =1
nplt = floor((tmax/0.01)/dt)
nmax = int(round(tmax/dt))


for t in range(1,nmax+1):
    for q in range(1, m):
        p=q-1
        r=q+1
        if q==1:
            p=m
        if q==m:
            r=1
        A=3/2*onsarray-dt*F[q]*onsarray
        cons=1/A
        B=2*u_hat[q]-0.5*uold_hat[q]-gamma*dt/2*(1j*kk+ll)*u_hat[r]-gamma*dt/2*(1j*kk-ll)*unew_hat[q]
        unew_hat=B*cons
        uold_hat=u_hat 
        u_hat=unew_hat
        u1 = real(fft.ifft(u_hat))
        uold=real(fft.ifft(uold_hat)) 

        


