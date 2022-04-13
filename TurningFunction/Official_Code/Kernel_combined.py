from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special
from compute_3D_FFT import *
from Kernel_d import *
from Kernel_o import *


def Kernel(N,m,d,a_or_r):
    temp_2D = zeros((N,N))
    temp_3D=[]
    temp2_3D=[]
    for i in range (0,N):
        array(temp_3D.append(temp_2D))
        array(temp2_3D.append(temp_2D))
    K_d=temp_2D
    K_o=temp_3D
    Kernel_func=array(temp2_3D)
    result=array(temp2_3D)
    K_d=array(Kernel_d(N,m,d))
    K_o=array(Kernel_o(N,a_or_r))
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0,N):
                kd=array(K_d)[j,k]
                ko=array(K_o)[i,j,k]
                temp=kd*ko
                Kernel_func[i,j,k]=temp
    result=Kernel_func
    return result


