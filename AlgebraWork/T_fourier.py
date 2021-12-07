from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA

def T_hat(N,m,u_hat,kenel_d_hat,kernel_o_hat,turning_hat):
    u3D_hat=[]
    for i in range(0,m):
        u2D_hat_temp=u_hat
        for k in range(0,N*N):
            u2D_hat_temp[k]=u2D_hat_temp[k]*kenel_d_hat
        for j in range(0,m):
            u2D_hat_temp[:,j]=u2D_hat_temp[:,j]*kernel_o_hat*w_hat[i,:]
        u3D_hat.append(u2D_hat_temp)
        for k in range(0,N*N):
            for j in range(0,m):
                u3D_hat[:,i,j]=u3D_hat[:,i,j]*w_hat[:,i]
        return u3D_hat
            
