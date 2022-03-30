from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special


def K_w_combine(N, Kernel,w):
    Kernel=array(Kernel)
    w=array(w)
    temp_2D=array(zeros((N,N)))
    temp_3D=[]
    for i in range (0,N):
        array(temp_3D.append(temp_2D))
    temp_4D=[]
    for i in range (0,N):
        array(temp_4D.append(temp_3D))
    temp_4D=array(temp_4D)
    K_w=temp_4D
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0,N):
                for l in range (0,N):
                    k_temp=Kernel[j,k,l]
                    w_temp=w[i,j,k,l]
                    k_w_temp=k_temp*w_temp
                    K_w[i,j,k,l]=k_w_temp
    return K_w
        