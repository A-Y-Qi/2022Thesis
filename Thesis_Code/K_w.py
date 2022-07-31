from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special
from matplotlib import *
from compute_2D_IFFT import *

'''
Kernels
'''

def Kernel_d(N,m,d):
    #define x and y
    x = (2*pi/N)*arange(-N/2.0,N/2.0)
    y = (2*pi/N)*arange(-N/2.0,N/2.0)
    #define constant A
    A=pi*m*(m*exp(-(d**2)/(m**2))+d*sqrt(pi)+d*sqrt(pi)*special.erf(d/m))
    #Generate a 2*2 matrix
    Kernel_d_func = zeros((N,N))
    for i in range (0,N):
        for j in range(0,N):
            Kernel_d_func[i,j]=(1/A)*exp((-(sqrt(x[j]**2+y[i]**2)-d)**2)/(m**2))
    return Kernel_d_func

def Kernel_o(N,a_or_r):
    #define theta and phi
    phi = (2*pi/N)*arange(-N/2.0,N/2.0)
    #define Sx and Sy
    pos = (2*pi/N)*arange(-N/2.0,N/2.0)   
    #Generate a 2*2 matrix
    Kernel_o_func=zeros((N,N,N))
    Kernel_o_func=array(Kernel_o_func)
    for i in range (0,N):
        for k in range(0,N):
            for l in range (0,N):
                psi=arctan2(pos[k],pos[l])
                temp_cos=cos(psi)
                temp_sin=sin(psi)
                if a_or_r=="a":
                    temp=(1/(2*pi))*(-(cos(phi[i])*temp_cos+sin(phi[i])*temp_sin)+1)
                elif a_or_r=="r":
                    temp=(1/(2*pi))*(cos(phi[i])*temp_cos+sin(phi[i])*temp_sin+1)
                Kernel_o_func[i,k,l]=temp
    return Kernel_o_func

def Kernel(N,a_or_r):
    result = zeros((N,N,N))
    #define theta and phi
    phi = (2*pi/N)*arange(-N/2.0,N/2.0)
    #define Sx and Sy
    space = (2*pi/N)*arange(-N/2.0,N/2.0)
    if a_or_r=="r":
        #my value
        #d=0
        #Lennaert value
        d=0
        m=0.1
        A=pi*m*(m*exp(-(d**2)/(m**2))+d*sqrt(pi)+d*sqrt(pi)*special.erf(d/m))
        kk=1
    elif a_or_r=="a":
        #my value
        #d=pi/4
        #Lennaert value
        d=0.5
        m=0.5
        A=pi*m*(m*exp(-(d**2)/(m**2))+d*sqrt(pi)+d*sqrt(pi)*special.erf(d/m))
        kk=-1
    for i in range (0,N):
        for k in range (0,N):
            for l in range (0,N):
                K_d=(1/A)*exp((-(sqrt(space[l]**2+space[k]**2)-d)**2)/(m**2))
                K_o=(1/(2*pi))*(kk*(cos(phi[i]-arctan2(space[k],space[l])))+1)
                result[i,k,l]=K_d*K_o
    return result




'''
periodic gaussian function
'''
# Approximation to the periodic Gaussian (truncated at 2*tru+1 terms):
def g(theta,sigma,N):
    pg = 0.0
    for z in range(-N,N+1):
        pg += exp(-((theta+2.0*pi*float(z))/sigma)**2)
    pg /= sigma*sqrt(pi)
    return pg

'''
periodic Gaussian probability function
'''
def w_probability_g(N, sigma,kk,a_or_r):
    #define theta and phi
    phi_prime = (2*pi/N)*arange(-N/2.0,N/2.0)
    phi = (2*pi/N)*arange(-N/2.0,N/2.0)
    #define Sx and Sy
    pos = (2*pi/N)*arange(-N/2.0,N/2.0)
    #a or r
    if a_or_r=="a":
        kk=abs(kk)
    if a_or_r=="r":
        kk=-abs(kk)
    #set up frame
    w_func_2D = array(zeros((N,N)))
    w_func_3D=array(zeros((N,N,N)))
    w_func_4D=array(zeros((N,N,N,N)))
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0,N):
                for l in range (0,N):
                    psi=arctan2(pos[k],pos[l])
                    v_angle=sin(phi_prime[j]-psi)
                    v=kk*v_angle 
                    phi_angle=phi_prime[j]-phi[i]                
                    theta=phi_angle-v
                    w_func_4D[i,j,k,l]=g(theta, sigma,N)
    return w_func_4D

    

#=====================Kw_combined function =============================#
#All variables are defined inside of the functions
#m,d are defined in Kernel
#sigma kk are defined in K_w


'''
K_w function with periodic Gaussian probability function
'''
def K_w_g(N,a_or_r):
    sigma=2.0
    kk=1.0
    Kernel_func=Kernel(N,a_or_r)
    #w_func=w_probability(N, sigma,kk,a_or_r)
    w_func=w_probability_g(N, sigma,kk,a_or_r)
    Kw_func_4D=copy(w_func)
    for i in range (0,N):
        for j in range (0,N):
            Kw_func_4D[i,j,:,:]=Kernel_func[j,:,:]*w_func[i,j,:,:]            
    return Kw_func_4D  

'''
K_w function with periodic Gaussian probability function in Fourier space
'''
def K_w_g_fft(N,K_w_g):
    Kw_func_4D_fft=zeros((N,N,N,N)).astype(complex)
    sign=sign_M(N)
    Kw_func_4D_fft=fft_2D_4DM(K_w_g,sign, N)
    return Kw_func_4D_fft

