# Time-stepping for "Fetecau with saturation" using the ETD1 method.
import numpy as np
from time import time                               # Just to check the wall time.
from computeInteractions import computeInteractions # Compute the interaction terms.
from r2f_f2r import *                               # Fourier->real->Fourier for space and angle.
from plotSolution import plotSolution               # Plot density and, if arrows=1, arrows on top.
import matplotlib # The next four lines are for the animation:
import matplotlib.pyplot as plt
matplotlib.use("Agg")
from matplotlib.animation import FFMpegWriter       # If youdo not have ffmpeg installed, you can change to, e.g., ImageMagickWriter to create GIFs instead.

def march(ps,pn,pni,u,A,B,Qa,Qr,Ql):
    # Extract some parameters to make the code more readable:
    h = pn[0]
    nsteps = pni[0]
    n = pni[1]
    m = pni[2]
    mass = ps[16]
    
    # Set up the animation:
    metadata = dict(title='Aggregation', artist='Matplotlib',comment='See plotSolution.py for details')
    writer = FFMpegWriter(fps=15, metadata=metadata)
    fig = plt.figure() # Open figure with two sets of axes for the density plot and the colorbar.
    ax1 = fig.add_axes([0.1,0.1,0.7,0.8])
    ax2 = fig.add_axes([0.85,0.1,0.075,0.8])


    # Initialize time and output array:
    t = 0.0
    out = np.empty((nsteps,m+2))
    print('Computing %d time steps of size %e...' % (nsteps,h))

    # Linear part is propagated on a mixed Fourier-grid basis:
    uh = r2f_space(n,m,u)
    wtime = time()
    with writer.saving(fig, "aggregation_test.mp4", 100): # Write frames to the animation on the fly!
        for i in range(nsteps): # Loop over time steps
            Nh, Nmin, Nmax = computeInteractions(n,m,uh,Qa,Qr,Ql,ps,pni) # Compute interaction terms.
            # Actual ETD1 step:
            uh = A*uh + B*Nh
            t += h
            # Compute the total mass for each angle - should be constant without interactions:
            out[i,0] = t
            for j in range(1,m+1):
                out[i,j] = (2.0*np.pi)**2 * np.real(uh[0,0,j-1])
            out[i,m+1] = np.sum(out[i,1:m+1]) # Total mass - should be constant over time to time-stepping accuracy.
            u = f2r_space(n,m,uh)
            print('Range of u is [%e,%e], range of T is [%e,%e].' % (np.min(u),np.max(u),Nmin,Nmax))
            plotSolution(n,m,u,mass,fig,ax1,ax2,pni)          # Draw a movie frame.
            writer.grab_frame()
            ax1.clear()
            ax2.clear()

            # Echo progress:
            if np.mod(i,np.min([50,int(nsteps/10)])) == 0:
                print('%d percent done!' % (int(100*i/nsteps)))
    wtime = time() - wtime
    print('Wall time for time stepping was %f.' % (wtime))
    File = open("last","w")
    np.savetxt(File,np.reshape(u,(n*n*m,1)),fmt='%16.9e',delimiter=' ')

    return out

