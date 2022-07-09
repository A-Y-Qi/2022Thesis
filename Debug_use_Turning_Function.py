from numpy import *
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from numpy import linalg as LA
from scipy import special
from compute_2D_FFT import *
from compute_3D_FFT import *
from compute_2D_IFFT import *
from Kw_combined_debug import *
from S_u_combine_debug import *
from scipy import signal as sig

#q: a constant
#K_w_fft: the K_w function with its fft form over X, ie K_w_fft=fft_2D_4DM(K_w_func,N)
#turning_fun: 2D for alignment and 4D for attraction/repulsion
#u: a 3D function

def Turning(q, sign_matrix, K_w_fft,u, N):
    ##create 4D array
    temp_2D=zeros((N,N))
    temp_3D=zeros((N,N,N))
    Fourier_turning=array(zeros((N,N,N,N))).astype(complex)
    turning=array(zeros((N,N,N,N)))
    ##Integral over theta
    sum_phi=copy(temp_2D).astype(complex)
    for k in range (0,N):
        for l in range (0,N):
            sum_phi[k,l]=(2*pi/N)*real(fft.fft(u[:,k,l])[0])
    sum_phi_fft=fft.fft2(sum_phi)*sign_matrix*((1/N)**2)
    #wavenumber for phi
    for i in range (0,N):
        #wavenumber for phi'
        for j in range (0,N):
            Fourier_turning[i,j,:,:]=K_w_fft[i,j,:,:]*sum_phi_fft*sign_matrix
            turning[i,j,:,:]=q*((N)**2)*real(fft.ifft2(Fourier_turning[i,j,:,:]))
    #print('range of T is [%e,%e].' % (amin(turning),amax(turning)))
    return turning




def Turning_3(q, K_w_fft,u, N):
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
    #K_w_func=K_w(N,a_or_r)
    #K_w_fft=fft_2D_4DM(K_w_func,N)
    fft_u=fft_3D(u,N)
    ##Integral over theta
    sum_phi=copy(temp_2D).astype(complex)
    for k in range (0,N):
        for l in range (0,N):
            sum_phi[k,l]=(2*pi/N)*real(fft.fft(u[:,k,l])[0])
    sum_phi_fft=fft.fft2(sum_phi)
    #wavenumber for phi
    for i in range (0,N):
        #wavenumber for phi'
        for j in range (0,N):
            Fourier_turning[i,j,:,:]=K_w_fft[i,j,:,:]*sum_phi_fft*sign_matrix
            turning[i,j,:,:]=fft.fftshift(real(fft.ifft2(Fourier_turning[i,j,:,:])))
            #wavenumber for x
            #for k in range (0, N):
                ##wavenumber for y
                #for l in range (0, N):
                    #Kw_fft_temp=K_w_fft[i,j,k,l]
                    #u_fft_temp=sum_phi[k,l]                    
                    ##temp2=q*a_kenel_fft[j,k,l]*turn_fft[i,plus,k,l]*fft_u[int(N/2),k,l]
                    ##testuse
                    #fft_temp=q*Kw_fft_temp*u_fft_temp
                    #Fourier_turning[i,j,k,l]=fft_temp
    #turning=ifft_2D_4DM(Fourier_turning,N)
    #for i in range (0,N):
        ##wavenumber for phi'
        #for j in range (0,N):
            #turning[i,j,:,:]=fft.fftshift(fft.ifft2(Fourier_turning[i,j,:,:]))
    #return turning,Fourier_turning
    return turning

def Turning_2(q,K_w, u, N):
   ##create 4D array
    temp_2D=zeros((N,N))
    temp_3D=[]
    for i in range (0,N):
        array(temp_3D.append(temp_2D))
    turning=[]
    for i in range (0,N):
        array(turning.append(temp_3D))
    turning=array(turning).astype("float64")
    phi_sum=copy(temp_2D)
    for k in range (0,N):
        for l in range (0,N):
            phi_sum[k,l]=(2*pi/N)*real(fft.fft(u[:,k,l])[0])
    for i in range (0,N):
        for j in range (0,N):
            turning[i,j,:,:]=((2*pi/N)**2)*sig.convolve2d(phi_sum,K_w[i,j,:,:],mode='same')
    return turning
'''
K_w function
'''
#N=32
#sigma=pi/8
#gamma=1
##Fkk=1
##test on K_attract
##test3=fft_3D(K_attract,N)
#position=(2*pi/N)*arange(-N/2.0,N/2.0)
#X,Y=meshgrid(position, position)
#dt=0.01
##testing on the w function
##test2=fft_4D(w,N)
#test_angle=0
#q=1
#sample_2D=zeros((N, N))
#sign_matrix=copy(sample_2D)
#w_number_sign=(-1)**(arange(-N/2.0,N/2.0)+1)
#for i in range (0,N):
    #sign_matrix[:,i]=w_number_sign[i]*w_number_sign
#K_w_a=load("K_w_attraction_09.npy")
#K_w_r=load("K_w_repulsion_09.npy")
#K_w_fft_a=fft_2D_4DM(K_w_a,sign_matrix, N)
#K_w_fft_r=fft_2D_4DM(K_w_r,sign_matrix, N)
#save("K_w_fft_a", K_w_fft_a)
#save("K_w_fft_r", K_w_fft_r)


'''
u set up
'''
##Construct initial condition
#sample_2D=zeros((N, N))
#sample_3D=[]
#gaussian=exp(-(position**2)/(pi/16))
#for i in range (0,N):
    #sample_3D.append(sample_2D)
#sample_3D=array(sample_3D) #Will use for construct u0
#sample_3D_2=copy(sample_3D)
#sample_3D_3=copy(sample_3D)
#Non_lin=sample_3D_2
#sample_4D=[]
#for i in range (0,N):
    #sample_4D.append(sample_3D)
#sample_4D=array(sample_4D).astype("float64")

###Construct a testing matrix
#h=2*pi/N
#i=16
#for i in range (0,N):
    #for k in range (0,N):
        #for l in range(0,N):
            #sample_3D[i,k,l]=1000*exp(-(position[k]**2)-(position[l]**2))
#u=array(sample_3D)

'''
Sum_phi
'''
#sum_phi=copy(sample_2D)
#for k in range (0,N):
    #for l in range (0,N):
        #sum_phi[k,l]=(2*pi/N)*real(fft.fft(u[:,k,l])[0])
##q=1
##u_linear_use=load("Fixed_fft_2D_3DM_debug_u.npy")
##u=u_linear_use[0]
##sum_phi=copy(sample_2D)
##for i in range (0,N):
    ##for k in range (0,N):
        ##for l in range (0,N):
            ##sum_phi[k,l]=real(fft.fft(u[:,k,l])[0])
    




#fig1 = plt.figure()
#ax1 = fig1.gca()
#angle=0
#angle_p=0
#con1=ax1.contourf(X, Y, sum_phi, cmap="coolwarm")
#fig1.colorbar(con1)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax1.set_xlabel('x')
#ax1.set_ylabel('y')
#plt.suptitle("Integral of u over phi function")



##==============plot use using convolution function ====================#
#turning_conv_a=Turning_2(q, K_w_a, u, N)
#turning_conv_r=Turning_2(q, K_w_r, u, N)
###turning_conv_r=load("turning_conv_r.npy")
#fig1 = plt.figure()
#ax1 = fig1.gca()
#angle=10
#angle_p=16
#con1=ax1.contourf(X, Y, turning_conv_r[angle,angle_p], cmap="coolwarm")
#fig1.colorbar(con1)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax1.set_xlabel('x')
#ax1.set_ylabel('y')
#plt.suptitle("The Repulsion Turning Function at angle position phi="+str(position[angle])[:6]+" phi'="+str(position[angle_p])[:6])


#turning_conv_a=load("turning_conv_a.npy")
#fig2 = plt.figure()
#ax2 = fig2.gca()
#angle=10
#angle_p=16
#con2=ax2.contourf(X, Y, turning_conv_a[angle,angle_p], cmap="coolwarm")
#fig2.colorbar(con2)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax2.set_xlabel('x')
#ax2.set_ylabel('y')
#plt.suptitle("The Attraction Turning Function at angle position phi="+str(position[angle])[:6]+" phi'="+str(position[angle_p])[:6])

'''
Convolution theorem plot and function
'''

#turning_convT_a=Turning(q, sign_matrix,K_w_fft_a,u, N)
#save("Test_use_turning_a", turning_convT_a)

#fig3 = plt.figure()
#ax3 = fig3.gca()
#angle=10
#angle_p=16
#con3=ax3.contourf(X, Y, turning_convT_a[angle,angle_p], cmap="coolwarm")
#fig3.colorbar(con3)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax3.set_xlabel('x')
#ax3.set_ylabel('y')
#plt.suptitle("The Attraction Turning Function at angle position phi="+str(position[angle])[:6]+" phi'="+str(position[angle_p])[:6])

#turning_convT_r=Turning(q, sign_matrix,K_w_fft_r,u, N)
#save("Test_use_turning_r", turning_convT_r)

#fig4 = plt.figure()
#ax4 = fig4.gca()
##angle=16
##angle_p=0
#con4=ax4.contourf(X, Y, turning_convT_r[angle,angle_p], cmap="coolwarm")
#fig4.colorbar(con4)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax4.set_xlabel('x')
#ax4.set_ylabel('y')
#plt.suptitle("The Repulsion Turning Function at angle position phi="+str(position[angle])[:6]+" phi'="+str(position[angle_p])[:6])





##==============plot use for sum phi ====================#

##Figure
#fig1 = plt.figure()
#ax1 = fig1.gca()
#con1=ax1.contourf(X, Y, u[0], cmap="coolwarm")
#fig1.colorbar(con1)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax1.set_xlabel('x')
#ax1.set_ylabel('y')
#plt.suptitle("The u function when phi=-pi")


#fig2 = plt.figure()
#ax2 = fig2.gca()
#con2=ax2.contourf(X, Y, test_Turning_r[angle,angle_p], cmap="coolwarm")
#fig2.colorbar(con2)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax2.set_xlabel('x')
#ax2.set_ylabel('y')
#plt.suptitle("The Saturation Repulsion Function at angle position phi="+str(position[angle])[:6]+" phi'="+str(position[angle_p])[:6])

#=====================Saturation Test use==========================#
#saturation=0.5*(1+tanh(test_Turning_a))
#fig3 = plt.figure()
#ax3 = fig3.gca()
#angle=16
#angle_p=20
#con3=ax3.contourf(X, Y, saturation[angle,angle_p], cmap="coolwarm")
#fig3.colorbar(con3)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax3.set_xlabel('x')
#ax3.set_ylabel('y')
#plt.suptitle("The Repulsion Turning Function at angle position phi="+str(position[angle])[:6]+" phi'="+str(position[angle_p])[:6])




#========================Animation Use ==============================#
#fig2 = plt.figure()
#ax2 = fig2.gca(projection='3d')
##ax2 = fig2.gca()
##timestep_number=99

#kk = 0
#def animate2(i):
    #global kk
    #Z2 = test_Turning[kk,0]
    #kk += 1
    #ax2.clear()
    #ax2.plot_surface(X,Y,Z2,rstride=1, cstride=1,cmap=plt.cm.jet,linewidth=0,antialiased=False)
    ##ax2.contourf(X,Y,Z2)
    #ax2.set_zlim(0,800)
    #ax2.set_xlabel('x')
    #ax2.set_ylabel('y')
    #plt.suptitle("Simulation using Matrix when N="+str(N)+" phi="+str(position[kk])+" dt="+str(dt))
    #if kk >= 31:
        #kk=0

#ax2.view_init(90, 90)
#anim2 = animation.FuncAnimation(fig2,animate2)
#plt.show()  


'''
Non_linear out 
'''
#Nout=Nonlin_out(N, turning_convT_a, u)
#Nin=Nonlin_in(N, turning_convT_a, u)
#fig4 = plt.figure()
#ax4 = fig4.gca()
##angle=16
##angle_p=0
#con4=ax4.contourf(X, Y, Nout[16], cmap="coolwarm")
#fig4.colorbar(con4)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax4.set_xlabel('x')
#ax4.set_ylabel('y')
#plt.suptitle("The Nonlinear out Function at angle position angle="+str(position[angle_p])[:6])
            
            
#fig5 = plt.figure()
#ax5 = fig5.gca()
##angle=16
##angle_p=0
#con5=ax5.contourf(X, Y, Nin[16], cmap="coolwarm")
#fig5.colorbar(con5)
##ax1.contourf(X,Y,Z)
##ax1.set_zlim(0,400)
#ax5.set_xlabel('x')
#ax5.set_ylabel('y')
#plt.suptitle("The Nonlinear in Function at angle position angle="+str(position[angle_p])[:6])

#fig5 = plt.figure()
#ax5 = fig5.gca()
#con5=ax5.contourf(X, Y, u[16], cmap="coolwarm")
#fig5.colorbar(con5)
#ax5.set_xlabel('x')
#ax5.set_ylabel('y')
#plt.suptitle("The Nonlinear in Function at angle position angle="+str(position[angle_p])[:6])

    
                    
            
         
                   
                    
                    
