from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special
from matplotlib import *

def saturation(T,N,ind):
    #when no saturation function applied
    if ind==0:
        saturation_T=copy(T)
    #these two functions are not applied in the results
    elif ind==1:
        H=5
        sig=0.1
        mu=5
        saturation_T=H*(pi/2.0 + arctan((T-mu)/sig))/pi
    elif ind==2:
        H=5
        sig=0.25
        mu=5
        saturation_T=H*(pi/2.0 + arctan((T-mu)/sig))/pi
    #The following functions are applied to the initial functions in results
    elif ind==3:
        H=2
        sig=1
        mu=2
        saturation_T=(1+tanh((T-mu)/sig))*(H/2)
    elif ind==4:
        H=5
        sig=1.2
        mu=5
        saturation_T=(1+tanh((T-mu)/sig))*(H/2)
    elif ind==5:
        H=6
        sig=1
        mu=6
        saturation_T=(1+tanh((T-mu)/sig))*(H/2)
    elif ind==6:
        H=6
        sig=2.5
        mu=6
        saturation_T=(1+tanh((T-mu)/sig))*(H/2)    
    return saturation_T



def Nonlin_in(N, S, u):
    temp_3D=array(zeros((N,N,N)))
    result=array(temp_3D)
    for i in range (0,N):
        for k in range (0,N):
            for l in range (0,N):
                result[i,k,l]=dot(S[i,:,k,l],u[:,k,l])*2.0*pi/float(N)
    return result

def Nonlin_out(N, S, u):
    temp_3D=array(zeros((N,N,N)))
    lamb=array(temp_3D)
    for j in range (0,N):
        for k in range (0,N):
            for l in range (0,N):
                lamb[j,k,l]=sum(S[:,j,k,l])*2.0*pi/float(N)*u[j,k,l]
    return -lamb

def Non_linear(N,T,u,ind):
    S=saturation(T,N,ind)
    result=u
    result=Nonlin_out(N, S, u)+Nonlin_in(N, S, u)
    return result