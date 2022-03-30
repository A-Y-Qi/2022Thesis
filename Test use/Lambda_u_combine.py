from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special


def Lambda_u_combine(N, Lambda, u):
    Lambda=array(Lambda)
    u=array(u)
    temp_2D=array(zeros((N,N)))
    temp_3D=[]
    for i in range (0,N):
        array(temp_3D.append(temp_2D))
    result=array(temp_3D)
    for j in range (0,N):
        for k in range (0,N):
            for l in range (0,N):
                temp_Lambda=Lambda[j,k,l]
                temp_u=u[j,k,l]
                result_temp=temp_Lambda*temp_u
                result[j,k,l]=result_temp
    return result
        
        