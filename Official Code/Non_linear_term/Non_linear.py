from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special
from matplotlib import *
from Turning_Function import *
from Kw_combined import *
from compute_2D_FFT import *

def saturation(T,N):
    #function centered at a
    a=50
    #Slope of the saturation function 
    b=1
    #highest value is 2*c
    c=50
    saturation_T=copy(T)
    saturation_T=(1+tanh((T-a)*b)))*c
    return saturation_T



def Nonlin_in(N, T, u):
    Saturation=array(saturation(T,N))
    u=array(u)
    temp_2D=array(zeros((N,N)))
    temp_3D=array(zeros((N,N,N)))
    temp_4D=array(zeros((N,N,N,N)))
    S_u=temp_4D
    for i in range (0,N):
        for j in range (0,N):
            S_u[i,j]=Saturation[i,j]*u[j]
    result=array(temp_3D)
    for i in range (0,N):
        for k in range (0,N):
            for l in range (0,N):
                result[i,k,l]=real((2*pi/N)*fft.fft(S_u[i,:,k,l])[0])
    return result

def Nonlin_out(N, T, u):
    Saturation=array(saturation(T,N))
    temp_2D=array(zeros((N,N)))
    temp_3D=array(zeros((N,N,N)))
    lamb=array(temp_3D)
    for i in range (0,N):
        for k in range (0,N):
            for l in range (0,N):
                lamb[i,k,l]=real((2*pi/N)*fft.fft(Saturation[:,i,k,l])[0])*u[i,k,l]
    return -lamb

def Non_linear(N,T,u):
    result=u
    result=Nonlin_out(N, T, u)+Nonlin_in(N, T, u)
    return result
