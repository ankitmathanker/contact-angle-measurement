import numpy as np
import matplotlib.pyplot as plt
from ase.io import read, write
from ase.visualize import view
import os

#function to read a frame and compute contact angle
def contact_angle(pdb):
    frame=read(pdb)
    #Find all O position
    Os=frame[frame.symbols=='O']
    #For each z-bin, find O at X or Y edges
    zmin=min(Os.positions[:,2])
    zmax=max(Os.positions[:,2])
    zbinedges=np.linspace(zmin, zmax, 15)  #numbe of bins used
    zbins=zbinedges[:-1]+(zbinedges[1]-zbinedges[0])/2
    Oz_xmin, Oz_xmax, Oz_ymin, Oz_ymax=[],[],[],[]
    for z in range(len(zbinedges)-1):
        Oz=Os[(Os.positions[:,2]<=zbinedges[z+1])*(Os.positions[:,2]>zbinedges[z])]
        Oz_xmin.append(min(Oz.positions[:,0]))
        Oz_xmax.append(max(Oz.positions[:,0]))
        Oz_ymin.append(min(Oz.positions[:,1]))
        Oz_ymax.append(max(Oz.positions[:,1]))
    #Interpolate and get a vector
    nfit=10   # used first 10 bins to approximate
    fitxmin=np.polyfit(Oz_xmin[:nfit], zbins[:nfit],deg=1)
    fitxmax=np.polyfit(Oz_xmax[:nfit], zbins[:nfit],deg=1)
    fitymin=np.polyfit(Oz_ymin[:nfit], zbins[:nfit],deg=1)
    fitymax=np.polyfit(Oz_ymax[:nfit], zbins[:nfit],deg=1)
    #Use vector to compute angle
    thetaxmin=np.arccos(np.dot([1, fitxmin[0]], [1,0])/np.linalg.norm([1, fitxmin[0]]))/np.pi*180
    thetaxmax=np.arccos(np.dot([-1, fitxmax[0]], [-1,0])/np.linalg.norm([-1, fitxmax[0]]))/np.pi*180
    thetaymin=np.arccos(np.dot([1, fitymin[0]], [1,0])/np.linalg.norm([1, fitymin[0]]))/np.pi*180
    thetaymax=np.arccos(np.dot([-1, fitymax[0]], [-1,0])/np.linalg.norm([-1, fitymax[0]]))/np.pi*180

    print(pdb + ' done')
    return np.array([thetaxmin, thetaxmax, thetaymin, thetaymax]), zbins, np.array([Oz_xmin, Oz_xmax, Oz_ymin, Oz_ymax]), np.array([fitxmin,fitxmax,fitymin,fitymax])

#Loop through frames
angles, zbins, Ozs, fits=[],[],[],[]
indframes=np.arange(100000, 1100000, 2500)
for i in indframes:
    angle, zbin, Oz, fit=contact_angle('frames_pure_water/frame_'+str(i)+'.pdb')
    angles.append(angle)
    zbins.append(zbin)
    Ozs.append(Oz)
    fits.append(fit)

#Plotting angle distribution and average angle
plt.rcParams["font.family"] = "Arial"
plt.rcParams["font.size"] = 20
fig, axes=plt.subplots(1,2,figsize=(9,4))

colors=['b','r','g','orange']
[axes[0].plot(range(len(indframes)),np.array(angles)[:,i],'o',c=colors[i]) for i in range(4)]
axes[0].axhline(y=np.mean(angles), c='k', ls='--')
#axes[0].set_title('Time series')
axes[0].tick_params(which='major', length=6, width=4, axis="y", direction="in")
axes[0].tick_params(which='major', length=6, width=4, axis="x", direction="in")
axes[0].spines[['right', 'left', 'bottom', 'top']].set_linewidth(4)
axes[0].xaxis.set_tick_params(pad=10)  # Adjust padding for x-axis
axes[0].yaxis.set_tick_params(pad=10)  # Adjust padding for y-axis

axes[1].hist(np.array(angles).flatten(),facecolor='dimgrey',bins=30)
axes[1].axvline(x=np.mean(angles), c='w', ls='--')
#axes[1].set_title('Histogram')
# [i.tick_params(axis='both',labelsize=14) for i in axes]

axes[1].tick_params(which='major', length=6, width=4, axis="y", direction="in")
axes[1].tick_params(which='major', length=6, width=4, axis="x", direction="in")
axes[1].spines[['right', 'left', 'bottom', 'top']].set_linewidth(4)
# axes[1].xaxis.set_tick_params(pad=10)  # Adjust padding for x-axis
# axes[1].yaxis.set_tick_params(pad=10)  # Adjust padding for y-axis

plt.savefig('frames_and_degree.svg')
plt.show()

print('Avg. contact angle: {:.2f} +- {:.2f}'.format(np.mean(angles), np.std(angles)/np.sqrt(len(angles))))