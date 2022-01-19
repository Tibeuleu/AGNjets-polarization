#!/usr/bin/python
import os
import matplotlib.offsetbox
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm

from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator
from mpl_toolkits.axes_grid1.anchored_artists import (AnchoredSizeBar)

##################

pas = 750
pas2 = 10

#dossier = "data/"
#dossier = "data_1cm/"
#dossier = "data_2cm/"
dossier = "data_6cm/"
#dossier = "data_20cm/"

A=np.loadtxt(dossier+"photons.res")

theta_det = A[0,0]
pas_x = A[0,1]
pas_y = A[0,2]
nu = A[0,3]
y_bulk_j = A[0,4]
theta_j = A[1,0]
R_j = A[1,1]
n_photons = A[1,2]
p_j = A[1,3]
B_j = A[1,4]
D_l = A[2,0]
z = A[2,1]
ne = A[2,2]
B_f = A[2,3]

for k in range(10):
    fig = plt.figure(figsize=(18,9))
    ax1 = fig.add_subplot(121,projection='3d')
    B=np.loadtxt(dossier+"photons_{}.res".format(k))

    X = A[3::pas,0]
    Y = A[3::pas,1]
    Z = A[3::pas,2]
    theta = B[::pas2,3]*np.pi/180
    phi = B[::pas2,4]*np.pi/180
    bx = np.sin(theta)*np.cos(phi)
    by = np.sin(theta)*np.sin(phi)
    bz = np.cos(theta)
    ax1.scatter(X,Y,Z)
    ax1.quiver(B[::pas2,0],B[::pas2,1],B[::pas2,2],bx,by,bz,color='r',length=2e20,normalize=True)
    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max() / 2.0
    mid_x = (X.max()+X.min()) * 0.5
    mid_y = (Y.max()+Y.min()) * 0.5
    mid_z = (Z.max()+Z.min()) * 0.5
    ax1.set_xlim(mid_x - max_range, mid_x + max_range)
    ax1.set_ylim(mid_y - max_range, mid_y + max_range)
    ax1.set_zlim(mid_z - max_range, mid_z + max_range)
    ax1.set_xlabel('Offset (m)')
    ax1.set_ylabel('Offset (m)')
    ax1.set_zlabel('Offset (m)')

    ax2 = fig.add_subplot(122)
    I=np.loadtxt(dossier+"detect_0_I_{}.res".format(k))

    im2 = ax2.imshow(np.asarray(np.log10(I[::-1,:])),extent=[-I.shape[1]/2., I.shape[1]/2., -I.shape[0]/2., I.shape[0]/2. ],vmin=0.,vmax=7.,aspect='auto',cmap='gnuplot2_r',alpha=0.2)
    divider2 = make_axes_locatable(ax2)
    cax2 = divider2.append_axes("right", size="2%", pad=0.05)
    cbar2 = plt.colorbar(im2, cax=cax2, label=r"$log_{10}$ Intensity (arbitrary unit)")
    cont2 = ax2.contour(np.log10(I),5,extent=[-I.shape[1]/2., I.shape[1]/2., -I.shape[0]/2., I.shape[0]/2. ],vmin=0.,vmax=7.,cmap='gnuplot2_r')
    ax2.clabel(cont2,inline=True,fontsize=8)
    ax2.set_xlabel('Offset (pixels)')
    ax2.set_ylabel('Offset (pixels)')
    ax2.axis("equal")

    ikpc = lambda x: x*3.085e19 #x in kpc, return in km
    bar = AnchoredSizeBar(ax2.transData, ikpc(10)/pas_x, '10 kpc', 3, pad=0.5, sep=5, borderpad=0.5, frameon=False, size_vertical=0.005, color='k')
    ax2.add_artist(bar)

    pola_p=np.loadtxt(dossier+"detect_0_pola_p_{}.res".format(k))
    pola_theta=np.loadtxt(dossier+"detect_0_pola_theta_{}.res".format(k))

    X, Y = np.linspace(-I.shape[0]/2.,I.shape[0]/2.,I.shape[0]/1.), np.linspace(-I.shape[1]/2.,I.shape[1]/2.,I.shape[1]/1.)
    xx, yy = np.meshgrid(X,Y)
    Q = ax2.quiver(xx[pola_p > 0],yy[pola_p > 0],(pola_p*np.cos(pola_theta))[pola_p > 0],(pola_p*np.sin(pola_theta))[pola_p > 0],pivot='mid',angles='xy', scale_units='xy', scale=1./1.,headwidth=0.)

    bar2 = AnchoredSizeBar(ax2.transData, 1., r"$p_{100\%}$", 4, pad=0.5, sep=5, borderpad=0.5, frameon=False, size_vertical=0.005, color='k')
    ax2.add_artist(bar2)

    plt.suptitle("Projected image of the jets for an ordered boosted synchrotron emission.\n"+"With phi = {}deg, jet aperture angle = {}deg, nu = {:0.0e} Hz, B = {:0.0e} T, gamma_max = 1e+4.5, p = {}, gamma_bulk = {}.".format(10*k,theta_j,nu,B_j,p_j,y_bulk_j))
    plt.savefig(dossier+"projection_{:02d}_0.png".format(int(k)),bbox_inches='tight')

os.system("convert -delay 50 -loop 0 "+dossier+"*_0.png "+dossier+"projection_log_0.gif")
#plt.show()
