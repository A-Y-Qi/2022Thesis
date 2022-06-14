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
    saturation_T=copy(T)
    #saturation_T=(1+tanh((T-700)/100))*0.5
    saturation_T=(1+tanh((T-amax(T)*0.5)*(1/(0.25*amax(T)))))*0.5
    return saturation_T



def S_u(N, T, u):
    Saturation=array(saturation(T,N))
    u=array(u)
    temp_2D=array(zeros((N,N)))
    temp_3D=[]
    for i in range (0,N):
        array(temp_3D.append(temp_2D))
    temp_4D=[]
    for i in range (0,N):
        array(temp_4D.append(temp_3D))
    temp_4D=array(temp_4D)
    S_u=temp_4D
    for i in range (0,N):
        for j in range (0,N):
            S_temp=Saturation[i,j]
            u_temp=u[j]
            S_u_temp=S_temp*u_temp
            S_u[i,j]=S_u_temp
    return S_u

def Nonlin_in(N, T, u):
    Saturation=array(saturation(T,N))
    u=array(u)
    temp_2D=array(zeros((N,N)))
    temp_3D=[]
    for i in range (0,N):
        temp_3D.append(temp_2D)
    temp_4D=[]
    for i in range (0,N):
        temp_4D.append(temp_3D)
    temp_4D=array(temp_4D)
    temp_3D=array(temp_3D)
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
    temp_3D=[]
    for i in range (0,N):
        temp_3D.append(temp_2D)
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
