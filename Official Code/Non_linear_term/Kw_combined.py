from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special
from matplotlib import *


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
            Kernel_d_func[i,j]=(1/A)*exp((-(sqrt(x[i]**2+y[j]**2)-d)**2)/(m**2))
    return Kernel_d_func

def Kernel_o(N,a_or_r):
    #define theta and phi
    phi = (2*pi/N)*arange(-N/2.0,N/2.0)
    #define Sx and Sy
    Sx = (2*pi/N)*arange(-N/2.0,N/2.0)
    Sy = (2*pi/N)*arange(-N/2.0,N/2.0)    
    #Generate a 2*2 matrix
    Kernel_o_func_1 = zeros((N,N))
    Kernel_o_func=[]
    for i in range (0,N):
        array(Kernel_o_func.append(Kernel_o_func_1))
    Kernel_o_func=array(Kernel_o_func)
    for i in range (0,N):
        for k in range(0,N):
            for l in range (0,N):
                if Sx[l]==Sy[k] and Sx[k]==0:
                    temp_sin=sin(phi[i])
                    temp_cos=cos(phi[i])
                else:
                        temp_cos=Sx[l]/sqrt(Sx[l]**2+Sy[k]**2)
                        temp_sin=Sy[k]/sqrt(Sx[l]**2+Sy[k]**2)
                if a_or_r=="a":
                    temp=(1/(2*pi))*(-cos(phi[i])*temp_cos-sin(phi[i])*temp_sin+1)
                elif a_or_r=="r":
                    temp=(1/(2*pi))*(cos(phi[i])*temp_cos+sin(phi[i])*temp_sin+1)
                Kernel_o_func[i,k,l]=temp
    return Kernel_o_func

def Kernel(N,a_or_r):
    temp_2D = zeros((N,N))
    temp_3D=[]
    for i in range (0,N):
        temp_3D.append(temp_2D)
    temp_3D=array(temp_3D)
    temp2_3D=copy(temp_3D)
    K_d=array(temp_2D)
    K_o=array(temp_3D)
    Kernel_func=array(temp2_3D)
    result=array(temp2_3D)
    m=pi/8
    if a_or_r=="r":
        d=0
    elif a_or_r=="a":
        d=pi/4
    K_d=array(Kernel_d(N,m,d))
    K_o=array(Kernel_o(N,a_or_r))
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0,N):
                kd=K_d[j,k]
                ko=K_o[i,j,k]
                temp=kd*ko
                Kernel_func[i,j,k]=temp
    result=Kernel_func
    return result
'''
Original Possibility function
(not correct)
'''

def w_possibility(N, sigma,kk,a_or_r):
    #define theta and phi
    phi_prime = (2*pi/N)*arange(-N/2.0,N/2.0)
    phi = (2*pi/N)*arange(-N/2.0,N/2.0)
    #define Sx and Sy
    Sx = (2*pi/N)*arange(-N/2.0,N/2.0)
    Sy = (2*pi/N)*arange(-N/2.0,N/2.0) 
    #a or r
    if a_or_r=="a":
        kk=abs(kk)
    if a_or_r=="r":
        kk=-abs(kk)
    #set up frame
    w_func_2D = array(zeros((N,N)))
    w_func_3D=[]
    for i in range (0,N):
        array(w_func_3D.append(w_func_2D))  
    w_func_4D=[]
    for i in range (0,N):
        array(w_func_4D.append(w_func_3D))
    w_func_4D=array(w_func_4D)
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0,N):
                for l in range (0,N):
                    Sxx=Sx[l]
                    Syy=Sy[k]
                    phi_prime_temp=phi_prime[j]
                    phi_temp=phi[i]
                    if Sxx==Syy and Sxx==0:
                        psi=phi_prime_temp   
                    else:
                        psi=arccos(Sxx/(sqrt((Sxx**2)+(Syy**2))))
                    if Syy<0:
                        psi=-psi
                    v_angle=phi_prime_temp-psi 
                    if v_angle>=pi:
                        v_angle=v_angle-2*pi
                    elif v_angle<=-pi:
                        v_angle=v_angle+2*pi
                    v=kk*v_angle 
                    phi_angle=phi_prime_temp-phi_temp
                    if phi_angle>=pi:
                        phi_angle=phi_angle-2*pi
                    elif phi_angle<=-pi:
                        phi_angle=phi_angle+2*pi                    
                    theta=phi_angle-v
                    if theta>=pi:
                        theta=theta-2*pi
                    elif theta<=-pi:
                        theta=theta+2*pi
                    if abs(theta)<=sigma:
                        w_func_4D[i,j,k,l]=1/(2*sigma)
                    else:
                        w_func_4D[i,j,k,l]=0
    return w_func_4D

'''
New possibility function
the one integrates to 1
'''

def w_possibility_2(N, sigma,kk,a_or_r):
    #define theta and phi
    phi_prime = (2*pi/N)*arange(-N/2.0,N/2.0)
    phi = (2*pi/N)*arange(-N/2.0,N/2.0)
    #define Sx and Sy
    Sx = (2*pi/N)*arange(-N/2.0,N/2.0)
    Sy = (2*pi/N)*arange(-N/2.0,N/2.0) 
    #a or r
    if a_or_r=="a":
        kk=abs(kk)
    if a_or_r=="r":
        kk=-abs(kk)
    #set up frame
    w_func_2D = array(zeros((N,N)))
    w_func_3D=[]
    for i in range (0,N):
        array(w_func_3D.append(w_func_2D))  
    w_func_4D=[]
    for i in range (0,N):
        array(w_func_4D.append(w_func_3D))
    w_func_4D=array(w_func_4D)
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0,N):
                for l in range (0,N):
                    Sxx=Sx[l]
                    Syy=Sy[k]
                    phi_prime_temp=phi_prime[j]
                    phi_temp=phi[i]
                    if Sxx==0 and Syy==0:
                        v_angle=phi_prime_temp
                    else:
                        v_angle=sin(phi_prime_temp)*(Sxx/sqrt(Sxx**2+Syy**2))-cos(phi_prime_temp)*(Syy/sqrt(Sxx**2+Syy**2))
                    v=kk*v_angle 
                    phi_angle=phi_prime_temp-phi_temp
                    if phi_angle>=pi:
                        phi_angle=phi_angle-2*pi
                    elif phi_angle<=-pi:
                        phi_angle=phi_angle+2*pi                    
                    theta=phi_angle-v
                    if theta>=pi:
                        theta=theta-2*pi
                    elif theta<=-pi:
                        theta=theta+2*pi
                    if abs(theta)<=sigma:
                        w_func_4D[i,j,k,l]=1/(2*sigma)
                    else:
                        w_func_4D[i,j,k,l]=0
    #for j in range (0,N):
        #for k in range (0,N):
            #for l in range (0,N):
                #if real(fft.fft(w_func_4D[:,j,k,l])[0])!=0:
                    #w_func_4D[:,j,k,l]=w_func_4D[:,j,k,l]*(1/real((2*pi/N)*(fft.fft(w_func_4D[:,j,k,l])[0])))
    return w_func_4D

def fft_2D_4DM(matrixin4D,sign_matrix, N):
        ##set up a new 3d matrix
        matrix_fft=array(matrixin4D).astype(complex)
        for i in range (0,N):
                for k in range (0,N):
                        #matrix_fft[i,k,:,:]=(((2*pi/N))**2)*fft.fftshift(fft.fft2(matrixin4D[i,k,:,:]))
                        matrix_fft[i,k,:,:]=(((2*pi/N))**2)*fft.fft2(matrixin4D[i,k,:,:])*sign_matrix
                
        return matrix_fft

    

#=====================Kw_combined function =============================#
#All variables are defined inside of the functions
#m,d are defined in Kernel
#sigma kk are defined in K_w

'''
Original K_w function
using original w function
'''

def K_w(N,a_or_r):
    sigma=pi/8
    kk=1
    Kernel_func=Kernel(N,a_or_r)
    #w_func=w_possibility(N, sigma,kk,a_or_r)
    w_func=w_possibility(N, sigma,kk,a_or_r)
    Kw_func_4D=copy(w_func)
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0,N):
                for l in range (0,N):
                    Kw_func_4D[i,j,k,l]=Kernel_func[j,k,l]*w_func[i,j,k,l]
    return Kw_func_4D  

'''
the K_w function
using new w function
'''

def K_w_2(N,a_or_r):
    sigma=pi/8
    kk=1
    Kernel_func=Kernel(N,a_or_r)
    #w_func=w_possibility(N, sigma,kk,a_or_r)
    w_func=w_possibility_2(N, sigma,kk,a_or_r)
    Kw_func_4D=copy(w_func)
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0,N):
                for l in range (0,N):
                    Kw_func_4D[i,j,k,l]=Kernel_func[j,k,l]*w_func[i,j,k,l]
    return Kw_func_4D  
