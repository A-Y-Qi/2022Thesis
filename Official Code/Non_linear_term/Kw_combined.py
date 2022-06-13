from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special
from matplotlib import *
from compute_2D_FFT import *

def Kernel_d(N,m,d):
    #define x and y
    x = (2*pi/N)*arange(-N/2.0,N/2.0)
    y = (2*pi/N)*arange(-N/2.0,N/2.0)
    #define constant A
    A=m*(m*exp(-(d**2)/(m**2))+d+d*special.erf(d/m))
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
                    if a_or_r=="r":
                        psi=phi[i]
                    elif a_or_r=="a":
                        psi=phi[i]
                else:
                    psi=arccos(Sx[l]/(sqrt((Sx[l]**2)+(Sy[k]**2))))
                if Sy[k]<0:
                    psi=-psi
                diff=phi[i]-psi
                temp=0
                if a_or_r=="r":
                    temp=(1/(2*pi))*(cos(phi[i]-psi)+1)
                if a_or_r=="a":
                    temp=(1/(2*pi))*(-cos(phi[i]-psi)+1)
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
                kd=array(K_d)[j,k]
                ko=array(K_o)[i,j,k]
                temp=kd*ko
                Kernel_func[i,j,k]=temp
    result=Kernel_func
    return result


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



#=====================Kw_combined function =============================#
#All variables are defined inside of the functions
#m,d are defined in Kernel
#sigma kk are defined in K_w

#def K_w(N,a_or_r):
    #sigma=pi/8
    #kk=1
    #Kernel_func=Kernel(N,a_or_r)
    #w_func=w_possibility(N, sigma,kk,a_or_r)
    #Kw_func_4D=copy(w_func)
    #for i in range (0,N):
        #for j in range (0,N):
            #for k in range (0,N):
                #for l in range (0,N):
                    #Kw_func_4D[i,j,k,l]=Kernel_func[j,k,l]*w_func[i,j,k,l]
    #return Kw_func_4D  

def K_w(N,a_or_r):
    sigma=pi/8
    kk=1
    Kernel_func=Kernel(N,a_or_r)
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
    Kw_func_4D=[]
    for i in range (0,N):
        array(Kw_func_4D.append(w_func_3D))
    Kw_func_4D=array(Kw_func_4D)
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
                        Kw_func_4D[i,j,k,l]=1/(2*sigma)*Kernel_func[j,k,l]
                    else:
                        Kw_func_4D[i,j,k,l]=0   
    return Kw_func_4D    
    
def K_w_fft(N,a_or_r):
    sigma=pi/8
    kk=1
    Kernel_func=Kernel(N,a_or_r)
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
    Kw_func_4D=[]
    for i in range (0,N):
        array(Kw_func_4D.append(w_func_3D))
    Kw_func_4D=array(Kw_func_4D)
    Kw_func_4D_fft=array(Kw_func_4D).astype(complex)
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
                        Kw_func_4D[i,j,k,l]=1/(2*sigma)*Kernel_func[j,k,l]
                    else:
                        Kw_func_4D[i,j,k,l]=0
    Kw_func_4D_fft=fft_2D_4DM(Kw_func_4D,N)*sign_M(N)
    return Kw_func_4D_fft
    
   
    

##grid_setup
#N=32
#position = (2*pi/N)*arange(-N/2.0,N/2.0)
#X,Y=meshgrid(position,position)

##function

#Kw_attr=K_w(N,"a")
#Kw_rep=K_w(N,"a")
##Kw_attr=load("Kw_attr.npy")
##Kw_rep=load("Kw_rep.npy")
##w_rep=load("w_rep.npy")
##w_attr=load("w_attr.npy")
##K_attr=Kernel(N,"a")
##K_rep=Kernel(N,"r")


##========================test 1================================#
#angle=16
#angle_prime=0
###Figure
#fig1 = plt.figure()
#ax1 = fig1.gca()
#con1=ax1.contourf(X, Y, Kw_rep[angle,angle_prime], cmap="coolwarm")
#fig1.colorbar(con1)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax1.set_xlabel('x')
#ax1.set_ylabel('y')
#plt.suptitle("The Repulsion K_w function at angle position phi'="+str(position[angle_prime])[:6]+" phi=" + str(position[angle])[:6])

#angle=0
#fig2 = plt.figure()
#ax2 = fig2.gca()
#con2=ax2.contourf(X, Y, Kw_rep[angle,angle_prime], cmap="coolwarm")
#fig2.colorbar(con2)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax2.set_xlabel('x')
#ax2.set_ylabel('y')
#plt.suptitle("The Repulsion K_w function at angle position phi'="+str(position[angle_prime])[:6]+" phi=" + str(position[angle])[:6])

#fig3 = plt.figure()
#ax3 = fig3.gca()
#con3=ax3.contourf(X, Y, K_rep[angle_prime], cmap="coolwarm")
#fig3.colorbar(con3)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax3.set_xlabel('x')
#ax3.set_ylabel('y')
#plt.suptitle("The Repulsion Kernel at angle position phi'="+str(position[angle_prime])[:6])

##========================test 2================================#
#angle=12
#angle_prime=12
###Figure
#fig4 = plt.figure()
#ax4 = fig4.gca()
#con4=ax4.contourf(X, Y, Kw_rep[angle,angle_prime], cmap="coolwarm")
#fig4.colorbar(con4)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax4.set_xlabel('x')
#ax4.set_ylabel('y')
#plt.suptitle("The Repulsion K_w function at angle position phi'="+str(position[angle_prime])[:6]+" phi=" + str(position[angle])[:6])

#angle=28
#fig5 = plt.figure()
#ax5 = fig5.gca()
#con5=ax5.contourf(X, Y, Kw_rep[angle,angle_prime], cmap="coolwarm")
#fig5.colorbar(con5)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax5.set_xlabel('x')
#ax5.set_ylabel('y')
#plt.suptitle("The Repulsion K_w function at angle position phi'="+str(position[angle_prime])[:6]+" phi=" + str(position[angle])[:6])

##fig6 = plt.figure()
##ax6 = fig6.gca()
##con6=ax6.contourf(X, Y, w_attr[angle,angle_prime], cmap="coolwarm")
##fig6.colorbar(con6)
###ax1.contourf(X,Y,Z)
###ax1.set_zlim(0,400)
##ax6.set_xlabel('x')
##ax6.set_ylabel('y')
##plt.suptitle("The Repulsion Kernel at angle position phi'="+str(position[angle_prime])[:6])

###========================test 3================================#
#angle=24
#angle_prime=24
###Figure
#fig7 = plt.figure()
#ax7 = fig7.gca()
#con7=ax7.contourf(X, Y, Kw_attr[angle,angle_prime], cmap="coolwarm")
#fig7.colorbar(con7)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax7.set_xlabel('x')
#ax7.set_ylabel('y')
#plt.suptitle("The Repulsion K_w function at angle position phi'="+str(position[angle_prime])[:6]+" phi=" + str(position[angle])[:6])

#angle=8
#fig8 = plt.figure()
#ax8 = fig8.gca()
#con8=ax8.contourf(X, Y, Kw_attr[angle,angle_prime], cmap="coolwarm")
#fig8.colorbar(con8)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax8.set_xlabel('x')
#ax8.set_ylabel('y')
#plt.suptitle("The Repulsion K_w function at angle position phi'="+str(position[angle_prime])[:6]+" phi=" + str(position[angle])[:6])

#fig2 = plt.figure()
#ax2 = fig2.gca()
#angle=4
#con2=ax2.contourf(X, Y, Kernel_rep[angle], cmap="coolwarm")
#fig2.colorbar(con2)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax2.set_xlabel('x')
#ax2.set_ylabel('y')
#plt.suptitle("The Repulsion Kernel of the individual at angle position phi="+str(position[angle])[:6])

#fig3 = plt.figure()
#ax3 = fig3.gca()
#angle=8
#con3=ax3.contourf(X, Y, Kernel_rep[angle], cmap="coolwarm")
#fig3.colorbar(con3)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax3.set_xlabel('x')
#ax3.set_ylabel('y')
#plt.suptitle("The Repulsion Kernel of the individual at angle position phi="+str(position[angle])[:6])

#fig4 = plt.figure()
#ax4 = fig4.gca()
#angle=12
#con4=ax4.contourf(X, Y, Kernel_rep[angle], cmap="coolwarm")
#fig4.colorbar(con4)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax4.set_xlabel('x')
#ax4.set_ylabel('y')
#plt.suptitle("The Repulsion Kernel of the individual at angle position phi="+str(position[angle])[:6])

#fig5 = plt.figure()
#ax5 = fig5.gca()
#angle=16
#con5=ax5.contourf(X, Y, Kernel_rep[angle], cmap="coolwarm")
#fig5.colorbar(con5)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax5.set_xlabel('x')
#ax5.set_ylabel('y')
#plt.suptitle("The Repulsion Kernel of the individual at angle position phi="+str(position[angle])[:6])

#fig6 = plt.figure()
#ax6 = fig6.gca()
#angle=20
#con6=ax6.contourf(X, Y, Kernel_rep[angle], cmap="coolwarm")
#fig6.colorbar(con6)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax6.set_xlabel('x')
#ax6.set_ylabel('y')
#plt.suptitle("The Repulsion Kernel of the individual at angle position phi="+str(position[angle])[:6])

#fig7 = plt.figure()
#ax7 = fig7.gca()
#angle=24
#con7=ax7.contourf(X, Y, Kernel_rep[angle], cmap="coolwarm")
#fig7.colorbar(con7)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax6.set_xlabel('x')
#ax6.set_ylabel('y')
#plt.suptitle("The Repulsion Kernel of the individual at angle position phi="+str(position[angle])[:6])


##Animation use
#fig10 = plt.figure()
##ax10 = fig10.gca(projection='3d')
#ax10 = fig10.gca()
##timestep_number=99

#count = 0
#angle_prime=22
#levels = linspace(0,1,100)
#def animate2(i):
    #global count
    #Z2=Kw_attr[count,angle_prime]
    #count += 1
    #ax10.clear()
    ##ax10.plot_surface(X,Y,Z2,rstride=1, cstride=1,cmap=plt.cm.jet,linewidth=0,antialiased=False)
    #ax10.contourf(X,Y,Z2,levels=levels,cmap=plt.cm.turbo)
    ##ax10.set_zlim(0,3)
    #ax10.set_xlabel('x')
    #ax10.set_ylabel('y')
    ##ax1.view_init(90, 90)
    ##plt.suptitle("Simulation using Matrix when N="+str(N)+" phi="+str(position[kk])+" dt="+str(dt))
    #if count >= N:
        #count=0

##ax10.view_init(90, 90)
#anim10 = animation.FuncAnimation(fig10,animate2)
#plt.show()    


