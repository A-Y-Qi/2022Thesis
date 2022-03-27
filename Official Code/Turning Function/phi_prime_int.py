from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special


def phi_prime_int(N, S_u):
    S_u=array(S_u)
    temp_2D=array(zeros((N,N)))
    temp_3D=[]
    for i in range (0,N):
        array(temp_3D.append(temp_2D))
    result=array(temp_3D)
    for i in range (0,N):
        for k in range (0,N):
            for l in range (0,N):
                temp_array=S_u[i,:,k,l]
                fft_S_u=(2*pi/N)*fft.fft(temp_array)
                fft_S_u=array(fft_S_u)
                integral=fft_S_u[0]
                result_int=real(integral)
                result[i,k,l]=result_int
    return result
        