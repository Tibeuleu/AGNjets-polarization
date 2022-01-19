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

#dossier = "data/"
#dossier = "data_1cm/"
#dossier = "data_2cm/"
#dossier = "data_6cm/"
dossier = "data_20cm/"

pas_pola = 5

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
    ax1 = fig.add_subplot(121)

    I=np.loadtxt(dossier+"detect_I_{}.res".format(k))

    pola_p_0=np.loadtxt(dossier+"detect_0_pola_p_{}.res".format(k))
    pola_theta_0=np.loadtxt(dossier+"detect_0_pola_theta_{}.res".format(k))

    pola_p_1=np.loadtxt(dossier+"detect_1_pola_p_{}.res".format(k))
    pola_theta_1=np.loadtxt(dossier+"detect_1_pola_theta_{}.res".format(k))

    for i in range(len(pola_p_0[:,0])):
        for j in range(len(pola_p_0[0,:])):
            if(i%pas_pola != 0 or j%pas_pola != 0):
                pola_p_0[i,j] = 0
    for i in range(len(pola_p_1[:,0])):
        for j in range(len(pola_p_1[0,:])):
            if(i%pas_pola != 0 or j%pas_pola != 0):
                pola_p_1[i,j] = 0

    p_max = np.amax(pola_p_0[pola_p_0 > 0])

    #print(p_max)

    im1 = ax1.imshow(np.asarray(np.log10(I[::-1,:])),extent=[-I.shape[1]/2., I.shape[1]/2., -I.shape[0]/2., I.shape[0]/2. ],vmin=0.,vmax=7.,aspect='auto',cmap='gnuplot2_r',alpha=0.2)
    divider1 = make_axes_locatable(ax1)
    cax1 = divider1.append_axes("right", size="2%", pad=0.05)
    cbar1 = plt.colorbar(im1, cax=cax1, label=r"$log_{10}$ Intensity (arbitrary unit)")
    cont1 = ax1.contour(np.log10(I),5,extent=[-I.shape[1]/2., I.shape[1]/2., -I.shape[0]/2., I.shape[0]/2. ],vmin=0.,vmax=7.,cmap='gnuplot2_r')
    ax1.clabel(cont1,inline=True,fontsize=8)
    ax1.set_xlabel('Offset (pixels)')
    ax1.set_ylabel('Offset (pixels)')
    ax1.set_title('Without Faraday effects')
    ax1.axis("equal")

    ikpc = lambda x: x*3.085e19 #x in kpc, return in km
    bar1 = AnchoredSizeBar(ax1.transData, ikpc(10)/pas_x, '10 kpc', 3, pad=0.5, sep=5, borderpad=0.5, frameon=False, size_vertical=0.005, color='k')
    ax1.add_artist(bar1)

    X, Y = np.linspace(-I.shape[0]/2.,I.shape[0]/2.,I.shape[0]/1.), np.linspace(-I.shape[1]/2.,I.shape[1]/2.,I.shape[1]/1.)
    xx, yy = np.meshgrid(X,Y)
    Q = ax1.quiver(xx[pola_p_0 > 0],yy[pola_p_0 > 0],(pola_p_0*np.cos(-pola_theta_0))[pola_p_0 > 0],(pola_p_0*np.sin(-pola_theta_0))[pola_p_0 > 0],pivot='mid',angles='xy', scale_units='xy', scale=p_max/(pas_pola),headwidth=0.,width=0.001)

    bar1 = AnchoredSizeBar(ax1.transData, pas_pola/p_max, r"$p_{max}$="+"{}".format(p_max), 4, pad=0.5, sep=5, borderpad=0.5, frameon=False, size_vertical=0.005, color='k')
    ax1.add_artist(bar1)


    ax2 = fig.add_subplot(122)

    im2 = ax2.imshow(np.asarray(np.log10(I[::-1,:])),extent=[-I.shape[1]/2., I.shape[1]/2., -I.shape[0]/2., I.shape[0]/2. ],vmin=0.,vmax=7.,aspect='auto',cmap='gnuplot2_r',alpha=0.2)
    divider2 = make_axes_locatable(ax2)
    cax2 = divider2.append_axes("right", size="2%", pad=0.05)
    cbar2 = plt.colorbar(im2, cax=cax2, label=r"$log_{10}$ Intensity (arbitrary unit)")
    cont2 = ax2.contour(np.log10(I),5,extent=[-I.shape[1]/2., I.shape[1]/2., -I.shape[0]/2., I.shape[0]/2. ],vmin=0.,vmax=7.,cmap='gnuplot2_r')
    ax2.clabel(cont2,inline=True,fontsize=8)
    ax2.set_xlabel('Offset (pixels)')
    ax2.set_ylabel('Offset (pixels)')
    ax2.set_title('With Faraday effects.')
    ax2.axis("equal")

    bar = AnchoredSizeBar(ax2.transData, ikpc(10)/pas_x, '10 kpc', 3, pad=0.5, sep=5, borderpad=0.5, frameon=False, size_vertical=0.005, color='k')
    ax2.add_artist(bar)

    Q = ax2.quiver(xx[pola_p_1 > 0],yy[pola_p_1 > 0],(pola_p_1*np.cos(-pola_theta_1))[pola_p_1 > 0],(pola_p_1*np.sin(-pola_theta_1))[pola_p_1 > 0],pivot='mid',angles='xy', scale_units='xy', scale=p_max/(pas_pola),headwidth=0.,width=0.001)

    bar2 = AnchoredSizeBar(ax2.transData, pas_pola/p_max, r"$p_{max}$="+"{}".format(p_max), 4, pad=0.5, sep=5, borderpad=0.5, frameon=False, size_vertical=0.005, color='k')
    ax2.add_artist(bar2)

    plt.suptitle("Projected image of the jets for an ordered boosted synchrotron emission taking into acount Faraday effects.\n"+"With phi = {}deg, jet aperture angle theta_j = {}deg, nu = {:0.0e} Hz, B = {:0.0e} T, gamma_max = 1e+4.5, p = {}, gamma_bulk = {}.".format(10*k,theta_j,nu,B_j,p_j,y_bulk_j)+"\n Physical parameters for Faraday effects : "+r"$D_{l}$="+"{:0.2e} ".format(D_l)+r"$pc$, $z$="+"{}, ".format(z)+r"$n_{e}$="+"{:0.2e} ".format(ne)+r"$m^{-3}$, $B_{intergalactique}$="+"{:0.1} ".format(B_f)+r"$\mu G$")
    plt.savefig(dossier+"projection_{:02d}.png".format(int(k)),bbox_inches='tight')

os.system("convert -delay 50 -loop 0 "+dossier+"*.png "+dossier+"comp_faraday.gif")
#plt.show()
