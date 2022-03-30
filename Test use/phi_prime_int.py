from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special
from compute_4D_FFT import *
from compute_3D_IFFT import *

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
                result_int=linalg.norm(integral)
                result[i,k,l]=result_int
    return result
        
def phi_prime_int_2(N, S_u):
    S_u=array(S_u)
    temp_2D=array(zeros((N,N)))
    temp_3D=[]
    for i in range (0,N):
        array(temp_3D.append(temp_2D))
    result=array(temp_3D)
    S_u_fft=fft_4D(S_u,N)
    Turning_FFT=S_u_fft[:,int(N/2),:,:]
    result=real(ifft_3D(Turning_FFT,N))
    #for i in range (0,N):
        #for k in range (0,N):
            #for l in range (0,N):
                #temp_array=S_u[i,:,k,l]
                #fft_S_u=fft_4D(N,S_u)
                #if linalg.norm(temp_array)>300:
                    #print(linalg.norm(temp_array))
                    #print(fft_S_u)
                #fft_S_u=array(fft_S_u)
                #integral=fft_S_u[0]
                #result_int=linalg.norm(integral)
                #result[i,k,l]=result_int
    return result