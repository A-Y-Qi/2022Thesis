from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special
from compute_4D_FFT import *
from compute_4D_IFFT import *

def w_func(N, sigma,kk,a_or_r):
    #define theta and phi
    phi_prime = (2*pi/N)*arange(-N/2.0,N/2.0)
    phi = (2*pi/N)*arange(-N/2.0,N/2.0)
    #define Sx and Sy
    Sx = (2*pi/N)*arange(-N/2.0,N/2.0)
    Sy = (2*pi/N)*arange(-N/2.0,N/2.0) 
    #a or r
    if a_or_r=="a":
        kk=abs(kk)
    if a_or_r=="r":
        kk=-abs(kk)
    #set up frame
    w_func_2D = array(zeros((N,N)))
    w_func_3D=[]
    for i in range (0,N):
        array(w_func_3D.append(w_func_2D))  
    w_func_4D=[]
    for i in range (0,N):
        array(w_func_4D.append(w_func_3D))
    w_func_4D=array(w_func_4D)
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0,N):
                for l in range (0,N):
                    Sxx=x[l]
                    Syy=Sy[k]
                    phi_prime_temp=phi_prime[j]
                    phi_temp=phi[i]
                    if Sxx==Syy and Sxx==0:
                        psi=0   
                    else:
                        psi=arccos(Sxx/(sqrt((Sxx**2)+(Syy**2))))
                    v_angle=phi_prime_temp-psi 
                    if v_angle>=pi:
                        v_angle=v_angle-2*pi
                    elif v_angle<=-pi:
                        v_angle=v_angle+2*pi
                    v=kk*v_angle 
                    phi_angle=phi_prime_temp-phi_temp
                    if phi_angle>=pi:
                        phi_angle=phi_angle-2*pi
                    elif phi_angle<=-pi:
                        phi_angle=phi_angle+2*pi                    
                    theta=phi_angle-v
                    if theta>=pi:
                        theta=theta-2*pi
                    elif theta<=-pi:
                        theta=theta+2*pi
                    if abs(theta)<=sigma:
                        w_func_4D[i,j,k,l]=1/(2*sigma)
                    else:
                        w_func_4D[i,j,k,l]=0
    return w_func_4D


def w_func_gaussian(N, sigma,kk,a_or_r):
    #define theta and phi
    phi_prime = (2*pi/N)*arange(-N/2.0,N/2.0)
    phi = (2*pi/N)*arange(-N/2.0,N/2.0)
    #define Sx and Sy
    Sx = (2*pi/N)*arange(-N/2.0,N/2.0)
    Sy = (2*pi/N)*arange(-N/2.0,N/2.0)
    #setupkk
    if a_or_r=="a":
        kk=abs(kk)
    if a_or_r=="r":
        kk=-abs(kk)    
    #set up frame
    w_func_2D = array(zeros((N,N)))
    w_func_3D=[]
    for i in range (0,N):
        array(w_func_3D.append(w_func_2D))  
    w_func_4D=[]
    for i in range (0,N):
        array(w_func_4D.append(w_func_3D))
    w_func_4D=array(w_func_4D)
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0,N):
                for l in range (0,N):
                    ssx=Sx[l]
                    ssy=Sy[k]
                    psi=0
                    theta=0
                    phi_prime_temp=phi_prime[j]
                    phi_temp=phi[i]
                    if ssx==ssy and ssx==0:
                        psi=0   
                    elif ssy>0:
                        psi=arccos(ssx/(sqrt((ssx**2)+(ssy**2))))
                    else:
                        psi=-arccos(ssx/(sqrt((ssx**2)+(ssy**2))))
                    diff=phi_prime_temp-psi
                    if diff>pi:
                        diff=-2*pi+diff
                    if diff<-pi:
                        diff=2*pi+diff
                    v=kk*(diff)
                    angle_temp=phi_prime_temp-v
                    if angle_temp>pi:
                        angle_temp=-2*pi+angle_temp
                    if angle_temp<-pi:
                        angle_temp=2*pi+angle_temp
                    theta=angle_temp-phi_temp
                    if theta>pi:
                        theta=-2*pi+theta
                    if theta<-pi:
                        theta=2*pi+theta
                    w_func_4D[i,j,k,l]=(1/(sqrt(pi)*sigma))*exp(-((theta)/sigma)**2)
    return w_func_4D



