from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special


def S_u_combine(N, S, u):
    Saturation=array(S)
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
            for k in range (0,N):
                for l in range (0,N):
                    S_temp=Saturation[i,j,k,l]
                    u_temp=u[j,k,l]
                    S_u_temp=S_temp*u_temp
                    S_u[i,j,k,l]=S_u_temp
    return S_u
        