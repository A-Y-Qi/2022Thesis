from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special


def Lambda(N, S):
    S=array(S)
    temp_2D=array(zeros((N,N)))
    temp_3D=[]
    for i in range (0,N):
        array(temp_3D.append(temp_2D))
    result=array(temp_3D)
    for j in range (0,N):
        for k in range (0,N):
            for l in range (0,N):
                temp_array=S_u[:,j,k,l]
                fft_S=(2*pi/N)*fft.fft(temp_array)
                fft_S=array(fft_S_u)
                integral=fft_S[0]
                result_int=real(integral)
                result[j,k,l]=result_int
    return result
        