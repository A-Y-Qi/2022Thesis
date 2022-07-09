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
            Kernel_d_func[i,j]=(1/A)*exp((-(sqrt(x[j]**2+y[i]**2)-d)**2)/(m**2))
    return Kernel_d_func

def Kernel_o(N,a_or_r):
    #define theta and phi
    phi = (2*pi/N)*arange(-N/2.0,N/2.0)
    #define Sx and Sy
    Sx = (2*pi/N)*arange(-N/2.0,N/2.0)
    Sy = (2*pi/N)*arange(-N/2.0,N/2.0)    
    #Generate a 2*2 matrix
    Kernel_o_func=zeros((N,N,N))
    Kernel_o_func=array(Kernel_o_func)
    for i in range (0,N):
        for k in range(0,N):
            for l in range (0,N):
                if l==k and k==N/2:
                    temp_sin=sin(phi[i])
                    temp_cos=cos(phi[i])
                else:
                        temp_cos=Sx[l]/sqrt(Sx[l]**2+Sy[k]**2)
                        temp_sin=Sy[k]/sqrt(Sx[l]**2+Sy[k]**2)
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
Wrong Version
'''
def w_probability(N, sigma,kk,a_or_r):
    #define theta and phi
    phi_prime = (2*pi/N)*arange(-N/2.0,N/2.0)
    phi = (2*pi/N)*arange(-N/2.0,N/2.0)
    #define Sx and Sy
    Sx = (2*pi/N)*arange(-N/2.0,N/2.0)
    Sy = (2*pi/N)*arange(-N/2.0,N/2.0) 
    #a or r
    if a_or_r=="a":
        kk=-abs(kk)
    if a_or_r=="r":
        kk=abs(kk)
    #set up frame
    w_func_2D = array(zeros((N,N)))
    w_func_3D=array(zeros((N,N,N)))
    w_func_4D=array(zeros((N,N,N,N)))
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
Correct version (delta function)
'''

def w_probability_2(N, sigma,kk,a_or_r):
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
    w_func_3D=array(zeros((N,N,N)))
    w_func_4D=array(zeros((N,N,N,N)))
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0,N):
                for l in range (0,N):
                    Sxx=Sx[l]
                    Syy=Sy[k]
                    phi_prime_temp=phi_prime[j]
                    phi_temp=phi[i]
                    if l==N/2 and k==N/2:
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
    return w_func_4D

'''
Correct version (periodic gaussian function)
'''
# Approximation to the periodic Gaussian (truncated at 2*tru+1 terms):
def g(theta,sigma,N):
    pg = 0.0
    for z in range(-N,N+1):
        pg += exp(-((theta+2.0*pi*float(z))/sigma)**2)
    pg /= sigma*sqrt(pi)
    return pg

def w_probability_g(N, sigma,kk,a_or_r):
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
    w_func_3D=array(zeros((N,N,N)))
    w_func_4D=array(zeros((N,N,N,N)))
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0,N):
                for l in range (0,N):
                    if l==N/2 and k==N/2:
                        v_angle=phi_prime[j]
                    else:
                        v_angle=sin(phi_prime[j])*(Sx[l]/sqrt(Sx[l]**2+Sy[k]**2))-cos(phi_prime[j])*(Sy[k]/sqrt(Sx[l]**2+Sy[k]**2))
                    v=kk*v_angle 
                    phi_angle=phi_prime[j]-phi[i]                
                    theta=phi_angle-v
                    w_func_4D[i,j,k,l]=g(theta, sigma,N)
    return w_func_4D

#fft 4D function
def sign_M(N):
    sample_2D=zeros((N, N))
    sign_matrix=copy(sample_2D)
    w_number_sign=fft.ifftshift((-1)**(arange(-N/2.0,N/2.0)+1))
    for i in range (0,N):
        sign_matrix[:,i]=w_number_sign[i]*w_number_sign
    return sign_matrix

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
Wrong version
'''
def K_w(N,a_or_r):
    sigma=pi/8
    kk=0.9
    Kernel_func=Kernel(N,a_or_r)
    w_func=w_probability(N, sigma,kk,a_or_r)
    Kw_func_4D=copy(w_func)
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0,N):
                for l in range (0,N):
                    Kw_func_4D[i,j,k,l]=Kernel_func[j,k,l]*w_func[i,j,k,l]
    return Kw_func_4D 

'''
Correct version
'''

def K_w_2(N,a_or_r):
    sigma=1.0
    kk=1.0
    Kernel_func=Kernel(N,a_or_r)
    #w_func=w_probability(N, sigma,kk,a_or_r)
    w_func=w_probability_2(N, sigma,kk,a_or_r)
    Kw_func_4D=copy(w_func)
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0,N):
                for l in range (0,N):
                    Kw_func_4D[i,j,k,l]=Kernel_func[j,k,l]*w_func[i,j,k,l]
    return Kw_func_4D  

'''
Periodic gaussian function
'''
def K_w_g(N,a_or_r):
    sigma=1.0
    kk=1.0
    Kernel_func=Kernel(N,a_or_r)
    #w_func=w_probability(N, sigma,kk,a_or_r)
    w_func=w_probability_g(N, sigma,kk,a_or_r)
    Kw_func_4D=copy(w_func)
    for i in range (0,N):
        for j in range (0,N):
            for k in range (0,N):
                for l in range (0,N):
                    Kw_func_4D[i,j,k,l]=Kernel_func[j,k,l]*w_func[i,j,k,l]
    return Kw_func_4D  

"""
Grid setup with N, position and XY
"""
N=64
#position = (2*pi/N)*arange(-N/2.0,N/2.0)
#X,Y=meshgrid(position,position)
#position_2 = (2/N)*arange(-N/2.0,N/2.0)
#X2,Y2=meshgrid(position_2,position_2)
#sigma=pi/8
#kk=1

#levels1=linspace(0,2.1,100)
#levels2=linspace(0,0.32,100)

'''
Distance Kernel
'''
#result_r=Kernel_d(N,pi/8,0)
#fig1 = plt.figure()
#ax1 = fig1.gca()
##ax1 = fig.gca()
##con8=ax8.contourf(X, Y, Z8,levels=levels,cmap=plt.cm.terrain)
#con1=ax1.contourf(X, Y, result_r,levels=levels1,cmap=plt.cm.terrain)
#fig1.colorbar(con1)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax1.set_xlabel('x')
#ax1.set_ylabel('y')
###ax6.view_init(90, 90)
#plt.suptitle("the Distance Kernel of Repulsion")

#result_a=Kernel_d(N,pi/8,pi/4)
#fig1 = plt.figure()
#ax1 = fig1.gca()
##ax1 = fig.gca()
##con8=ax8.contourf(X, Y, Z8,levels=levels,cmap=plt.cm.terrain)
#con1=ax1.contourf(X, Y, result_a,levels=levels2,cmap=plt.cm.terrain)
#fig1.colorbar(con1)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax1.set_xlabel('x')
#ax1.set_ylabel('y')
###ax6.view_init(90, 90)
#plt.suptitle("the Distance Kernel of Attraction")


'''
Produce and save the K_w function using K_w
'''
#K_w_a=K_w(N,"a")
#K_w_r=K_w(N,"r")
#K_w_a=load("K_w_64_attraction.npy")
#K_w_r=load("K_w_64_repulsion.npy")
#save("K_w_wrong_attraction_09",K_w_a)
#save("K_w_wrong_repulsion_09",K_w_r)

'''
Produce and save the K_w function using K_w_2
'''
#K_w_a=K_w_2(N,"a")
#K_w_r=K_w_2(N,"r")
#K_w_a=load("K_w_attraction.npy")
#K_w_r=load("K_w_repulsion.npy")
#save("K_w_64_attraction_Lennaert",K_w_a)
#save("K_w_64_repulsion_Lennaert",K_w_r)

'''
Produce and save the K_w function using K_w_g
'''
#sign_matrix=sign_M(64)
#K_w_a=K_w_2(N,"a")
#K_w_r=K_w_2(N,"r")
##K_w_a=load("K_w_attraction.npy")
##K_w_r=load("K_w_repulsion.npy")
#K_w_fft_a=fft_2D_4DM(K_w_a,sign_matrix, 64)
#K_w_fft_r=fft_2D_4DM(K_w_r,sign_matrix, 64)
#save("K_w_fft_64_attraction_Lennaert",K_w_fft_a)
#save("K_w_fft_64_repulsion_Lennaert",K_w_fft_r)


'''
plot K_w function 
'''
#angle=10
#angle_p=16
#fig1 = plt.figure()
#ax1 = fig1.gca()
##ax1 = fig.gca()
#Z1 = K_w_a[angle,angle_p]
##con8=ax8.contourf(X, Y, Z8,levels=levels,cmap=plt.cm.terrain)
#con1=ax1.contourf(X, Y, Z1,cmap=plt.cm.terrain)
#fig1.colorbar(con1)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax1.set_xlabel('x')
#ax1.set_ylabel('y')
###ax6.view_init(90, 90)
#plt.suptitle("the K_w function function with only attraction in effect at phi="+str(angle)+", phi_p="+str(angle_p))


#fig2 = plt.figure()
#ax2 = fig2.gca()
##ax1 = fig.gca()
#Z2 = K_w_r[angle,angle_p]
##con8=ax8.contourf(X, Y, Z8,levels=levels,cmap=plt.cm.terrain)
#con2=ax2.contourf(X, Y, Z2,cmap=plt.cm.terrain)
#fig2.colorbar(con2)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax2.set_xlabel('x')
#ax2.set_ylabel('y')
###ax6.view_init(90, 90)
#plt.suptitle("the K_w function function with only repulsion in effect at phi="+str(angle)+", phi_p="+str(angle_p))


'''
Plot probability function 
'''
#sigma=pi/8
#kk=1
##w_a=w_probability_2(N, sigma,kk,"a")
##w_r=w_probability_2(N, sigma,kk,"r")
##save("w_attraction",w_a)
##save("w_repulsion",w_r)
#w_a=load("w_attraction.npy")
#w_r=load("w_repulsion.npy")
#angle=0
#angle_p=16
#fig1 = plt.figure()
#ax1 = fig1.gca()
##ax1 = fig.gca()
#Z1 = w_a[angle,angle_p]
##con8=ax8.contourf(X, Y, Z8,levels=levels,cmap=plt.cm.terrain)
#con1=ax1.contourf(X, Y, Z1,cmap=plt.cm.terrain)
#fig1.colorbar(con1)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax1.set_xlabel('x')
#ax1.set_ylabel('y')
###ax6.view_init(90, 90)
#plt.suptitle("the w function function with only attraction in effect at phi="+str(angle)+", phi_p="+str(angle_p))


#fig2 = plt.figure()
#ax2 = fig2.gca()
##ax1 = fig.gca()
#Z2 = w_r[angle,angle_p]
##con8=ax8.contourf(X, Y, Z8,levels=levels,cmap=plt.cm.terrain)
#con2=ax2.contourf(X, Y, Z2,cmap=plt.cm.terrain)
#fig2.colorbar(con2)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax2.set_xlabel('x')
#ax2.set_ylabel('y')
###ax6.view_init(90, 90)
#plt.suptitle("the w function function with only repulsion in effect at phi="+str(angle)+", phi_p="+str(angle_p))





