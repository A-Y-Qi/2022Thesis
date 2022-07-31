from numpy import *

##In this function, we generate a 2D coefficient matrix
## with known x and y (l and k)
## where [I-dtM]^(-1)

##Make sure k and l are integers but not values
def diagonal_matrixB(N,k,l,gamma,dt):
    w_number= concatenate((linspace(0,N/2,int(N/2+1)),linspace(-N/2+1,-1,int(N/2-1))))
    w_number[int(N/2)] = 0
    kk=w_number[k]
    ll=w_number[l]     
    phi=(2*pi/N)*arange(-N/2.0,N/2.0)
    identity_matrix=identity(N)
    result=zeros(N).astype(complex)
    for i in range (0,N):
        theta=phi[i]
        M=-1*(ll*(cos(theta))+kk*(sin(theta)))
        const=dt*exp(1j*dt*gamma*M*0.5)*sinc(dt*0.5*M*gamma/pi)##note sinc in python has sin(pi*x)/(pi*x)
        result[i]=const
    return result

def matrixB(N,gamma,dt):
    result=array(zeros((N,N,N))).astype(complex)
    w_number= concatenate((linspace(0,N/2,int(N/2+1)),linspace(-N/2+1,-1,int(N/2-1))))
    w_number[int(N/2)] = 0
    phi=(2*pi/N)*arange(-N/2.0,N/2.0)
    for k in range (0,N):
        for l in range (0,N):
            M = - (w_number[l] * cos(phi) + w_number[k] * sin(phi))
            result[:,k,l]=dt * exp(0.5 * 1J * M * gamma * dt) * sinc(0.5 * M * gamma * dt / pi)
    return result
                
    
                
    
    
    
    