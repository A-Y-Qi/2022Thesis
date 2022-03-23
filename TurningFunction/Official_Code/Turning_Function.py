from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special
from compute_3D_FFT import *
from compute_2D_FFT import *
from compute_4D_FFT import *
from compute_3D_IFFT import *
from compute_2D_IFFT import *
from compute_4D_IFFT import *

#q: a constant
#dist_kernel: distance kernel, a 2D array
#orent_kernel: orientation kernel, a 1D array for alignment, 3D array for repulsion and attraction
#turning_fun: 2D for alignment and 4D for attraction/repulsion
#u: a 3D function
def Turning(q,K_w, u, N):
    ##create 4D array
    temp_2D=zeros((N,N))
    temp_3D=[]
    for i in range (0,N):
        array(temp_3D.append(temp_2D))
    Fourier_turning=[]
    turning=[]
    for i in range (0,N):
        array(Fourier_turning.append(temp_3D))
        array(turning.append(temp_3D))
    Fourier_turning=array(Fourier_turning).astype(complex)
    turning=array(turning).astype("float64")
    #Apply FFT to both equations K_w and u
    K_w_fft=fft_4D(K_w,N)
    fft_u=fft_3D(u,N)
    ##Integral over theta
    wave_0_u=fft_u[int(N/2)]
    for i in range (0,N):
        fft_u[i]=wave_0_u
    #wavenumber for phi
    for i in range (0,N):
        #wavenumber for phi'
        for j in range (0,N):
            #wavenumber for x
            for k in range (0, N):
                #wavenumber for y
                for l in range (0, N):
                    plus=mod(i+j,N)
                    Kw_fft_temp=K_w_fft[i,plus,k,l]
                    u_fft_temp=wave_0_u[k,l]                    
                    #temp2=q*a_kenel_fft[j,k,l]*turn_fft[i,plus,k,l]*fft_u[int(N/2),k,l]
                    #testuse
                    fft_temp=q*Kw_fft_temp*u_fft_temp
                    Fourier_turning[i,j,k,l]=fft_temp
    turning=ifft_4D(Fourier_turning,N)
    return turning

def Turning_2(q,K_w, u, N):
    ##create 4D array
    sample_2D=array(zeros((N,N)))
    complex_2D=array(zeros((N,N))).astype(complex)
    temp_3D=[]
    complex_3D=[]
    fft_u_3D=[]   
    for i in range (0,N):
        array(temp_3D.append(sample_2D))
        array(complex_3D.append(complex_2D))
        array(fft_u_3D.append(complex_2D)) 
    fft_u_3D=array(fft_u_3D)
    complex_4D=[]
    fft_Kw_4D=[]
    fft_turning=[]
    turning=[]
    for i in range (0,N):
        array(complex_4D.append(complex_3D))
        array(fft_Kw_4D.append(complex_3D))
        array(fft_turning.append(complex_3D))
        array(turning.append(temp_3D))        
    complex_4D=array(complex_4D)
    fft_Kw_4D=array(fft_Kw_4D).astype(complex)
    fft_turning=array(fft_turning).astype(complex)
    print(fft_turning[0,0,0,0])
    turning=array(turning) 
    #Apply 3D FFT to u
    fft_u=fft_3D(u,N)
    ##Integral over theta
    wave_0_u=fft_u[int(N/2)]
    for i in range (0,N):
        fft_u_3D[i]=wave_0_u
    fft_u=fft_u_3D
    #first calculate 2D fft for K_w
    K_w_2D=sample_2D
    K_w_fft_2D=complex_2D
    #wavenumber for phi
    for i in range (0,N):
        #wavenumber for phi'
        for j in range (0,N):
            K_w_2D=K_w[i,j]
            K_w_fft_2D=fft_2D(K_w_2D,N)
            fft_Kw_4D[i,j]=K_w_fft_2D
    #calculate fft_turnning
    #wavenumber for phi
    for i in range (0,N):
        #wavenumber for phi'
        for j in range (0,N):
            #wavenumber for x
            for k in range (0, N):
                #wavenumber for y
                for l in range (0, N):
                    Kw_fft_temp=fft_Kw_4D[i,j,k,l]
                    u_fft_temp=fft_u[j,k,l]
                    fft_temp=q*Kw_fft_temp*u_fft_temp
                    fft_turning[i,j,k,l]=fft_temp
                    fft_temp=0
    ##Find the 2D ifft for turning function
    #wavenumber for phi
    for i in range (0,N):
        #wavenumber for phi'
        for j in range (0,N):
            temp_2D=fft_turning[i,j]
            temp_ifft_2D=ifft_2D(temp_2D,N)
            turning[i,j]=temp_ifft_2D
    return turning
            
    
                    
            
         
                   
                    
                    
