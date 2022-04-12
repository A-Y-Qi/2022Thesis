from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from compute_3D_FFT import *
from compute_2D_FFT import *
from compute_3D_IFFT import *
from compute_2D_IFFT import *

def Moving(N,gamma,u):
    w_number=arange(-N/2.0,N/2.0)
    phi=(2*pi/N)*arange(-N/2.0,N/2.0)
    sample_2D=zeros((N,N))
    sample_3D=[]
    for i in range (0,N):
        sample_3D.append(sample_2D)
    sample_3D=array(sample_3D)
    sample_3D_2=array(sample_3D)
    sample_3D_3=array(sample_3D)
    sample_3D_4=array(sample_3D)
    u=array(u)
    u_x1=sample_3D_2
    u_x2=sample_3D_3
    result=sample_3D_4
    #find the cos(phi)*dx1_u term
    for i in range (0,N):
        for l in range (0,N):
            angle=phi[i]
            temp_x1=u[i,:,l]
            cos_phi=cos(angle)
            temp_fft_x1=fft.fftshift(fft.fft(temp_x1))
            fft_diff_x1=1j*w_number*temp_fft_x1
            ifft_diff_x1=real(fft.ifft(fft.ifftshift(fft_diff_x1)))
            u_x1[i,:,l]=ifft_diff_x1*cos_phi
    for i in range (0,N):
        for k in range (0,N):
            angle=phi[i]
            temp_x2=u[i,k,:]
            sin_phi=sin(angle)
            temp_fft_x2=fft.fftshift(fft.fft(temp_x2))
            fft_diff_x2=1j*w_number*temp_fft_x2
            ifft_diff_x2=real(fft.ifft(fft.ifftshift(fft_diff_x2)))
            u_x2[i,k,:]=ifft_diff_x2*sin_phi
    u_x1=array(u_x1)
    u_x2=array(u_x2)
    for i in range (0,N):
        for k in range (0,N):
            for l in range (0,N):
                x1_part=u_x1[i,k,l]
                x2_part=u_x2[i,k,l]
                result[i,k,l]=real(gamma*(x1_part+x2_part))
    return result
                
    
                
    
    
    
    
