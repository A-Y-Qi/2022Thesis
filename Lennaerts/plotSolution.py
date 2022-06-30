# Draw one frame of the movie. Total density at (x1,x2) with, if requested by setting arrows=1, arrows indicating local maxima over phi.
import numpy as np
import matplotlib.pyplot as plt
from auxiliary import setGrid

def extract_max(m,v,phis,thr):
    peaks = np.empty((m,2))
    nPeaks = 0
    vb = np.sum(v)/float(m)
    vs = np.sqrt(np.sum((v-vb)**2)/float(m))
    if v[m-1] < v[0] and v[0] > v[1] and v[0] > thr:
        peaks[nPeaks,0] = phis[0]
        peaks[nPeaks,1] = (v[0]-vb)/vs
        nPeaks += 1
    for j in range(1,m-1):
        if v[j-1] < v[j] and v[j] > v[j+1] and v[j] > thr:
            peaks[nPeaks,0] = phis[j]
            peaks[nPeaks,1] = (v[j]-vb)/vs
            nPeaks += 1
    if v[m-2] < v[m-1] and v[m-1] > v[0] and v[m-1] > thr:
        peaks[nPeaks,0] = phis[m-1]
        peaks[nPeaks,1] = (v[m-1]-vb)/vs
        nPeaks += 1
    peaks = peaks[0:nPeaks,:]
    return peaks
    
def plotSolution(n,m,u,mass,fig,ax1,ax2,pni):
    xs, phis = setGrid(pni)
    dum = np.sum(u,2) * 2.0*np.pi/float(m) # Integrate over the angle.
    # Note, that pcolormesh transposes the array by default.
    plotted = ax1.pcolormesh(np.transpose(dum))
    plt.colorbar(mappable=plotted, cax=ax2)

    if pni[4]:
        plt.sca(ax1)
        ub = np.sum(u)/float(n*n*m)
        us = np.sqrt(np.sum((u-ub)**2)/(n*n*m))
        thr = ub + 1.5*us
        skip = 4
        for i1 in range(0,n,skip):
            for i2 in range(0,n,skip):
                peaks = extract_max(m,u[i1,i2,:],phis,thr)
                for j in range(0,np.shape(peaks)[0]):
                    plt.arrow(i1,i2,peaks[j,1]*np.cos(peaks[j,0]),peaks[j,1]*np.sin(peaks[j,0]),color='white',width=0.1,head_width=0.5) 
    return
