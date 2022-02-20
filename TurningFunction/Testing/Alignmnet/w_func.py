from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special

def w_func(N, sigma,k):
    #define theta and phi
    theta = (2*pi/N)*arange(-N/2.0,N/2.0)
    phi = (2*pi/N)*arange(-N/2.0,N/2.0)
    v = k*theta
    pos_function=zeros((N,N))
    for i in range (0,N):
        for j in range (0,N):
            angle=phi[i]-v[j]
            if abs(angle)<=sigma:
                pos_function[i][j]=1/(2*sigma)
    return pos_function
    