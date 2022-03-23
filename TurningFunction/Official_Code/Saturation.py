from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special

def saturation_3D(T,N):
    sample_2D=array(zeros((N,N)))
    temp_3D=[]
    for i in range (0,N):
        array(temp_3D.append(sample_2D))
    saturation_T=array(temp_3D)
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0, N):   
                temp=0
                T_value=T[i,j,k]
                temp=0.5*(1+tanh(T_value))
                saturation_T[i,j,k]=temp
    return saturation_T

def saturation_4D(T,N):
    sample_2D=array(zeros((N,N)))
    sample_3D=[]
    for i in range (0,N):
        array(temp_3D.append(sample_2D))
    temp_4D=[]
    for i in range (0,N):
        array(temp_4D.append(sample_3D))    
    saturation_T=array(temp_4D)
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0, N):   
                for l in range (0,N):
                    temp=0
                    T_value=T[i,j,k,l]
                    temp=0.5*(1+tanh(T_value))
                    saturation_T[i,j,k,l]=temp
    return saturation_T
