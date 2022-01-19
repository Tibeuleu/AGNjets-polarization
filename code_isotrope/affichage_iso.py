#!/usr/bin/python
#-*- coding:utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

A=np.loadtxt("data_iso.res")

xdata_v = A[:,0]
xdata_l = 3e8/xdata_v
ydata_vi = A[:,1]
ydata_li = ydata_vi*3e8/(xdata_l**2)
ydata_vp = A[:,2]
ydata_lp = ydata_vp*3e8/(xdata_l**2)

fig,ax = plt.subplots(3,2,sharex='col',gridspec_kw={'height_ratios':[3,3,0.20]},figsize=(18,9))
plt.subplots_adjust(left=0.06,bottom=0.07,right=0.98,top=0.9,wspace=0.2,hspace=0.3)

ax[0,0].plot(xdata_l,ydata_li,label="for B = 10^-4 T,\n gamma = 10^4,\n pitch = 45deg.")
ax[0,0].set_xscale('log')
ax[0,0].set_xlabel(r"$\lambda$ ($m$)")
ax[0,0].set_yscale('log')
ax[0,0].set_xlim(3e-9,3e-1)
ax[0,0].set_ylabel(r"Emissivity ($W \cdot m^{-1}$)")
#ax[0,0].legend()
ax[0,0].set_title("Emissivity for one particle")

ax[0,1].plot(xdata_v,ydata_vi)
ax[0,1].set_xscale('log')
ax[0,1].set_yscale('log')
ax[0,1].set_xlim(1e9,1e17)
#ax[0,1].set_ylim(1e-8,1)
ax[0,1].set_xlabel(r"$\nu$ ($Hz$)")
ax[0,1].set_ylabel(r"Emissivity ($W \cdot Hz^{-1}$)")
#ax[0,1].legend()
ax[0,1].set_title("Emissivity for one particle")

ax[1,0].plot(xdata_l,ydata_lp)
ax[1,0].set_xlabel(r"$\lambda$ ($m$)")
ax[1,0].set_xscale('log')
ax[1,0].set_yscale('log')
ax[1,0].set_ylabel(r"Emissivity ($W \cdot m^{-1}$)")
#ax[1,0].legend()
ax[1,0].set_title("Emissivity for a population of particles")

ax[1,1].plot(xdata_v,ydata_vp)
ax[1,1].set_xscale('log')
ax[1,1].set_yscale('log')
#ax[1,1].set_ylim(1e9,1e17)
#ax[1,1].set_xlim(1e-2,1e2)
ax[1,1].set_xlabel(r"$\nu$ ($Hz$)")
ax[1,1].set_ylabel(r"Emissivity ($W \cdot Hz^{-1}$)")
#ax[1,1].legend()
ax[1,1].set_title("Emissivity for a population of particles")

def annotation_line(ax, xmin, xmax, y, text, ytext=0, linecolor='black', linewidth=1, fontsize=8):
    ax.annotate('', xy=(xmin, y), xytext=(xmax, y), xycoords='data', textcoords='data',arrowprops={'arrowstyle': '|-|', 'color':linecolor, 'linewidth':linewidth})
    xcenter = np.log10(xmin)+(np.log10(xmax)-np.log10(xmin))/2
    if ytext==0:
        ytext = y+(ax.get_ylim()[1]-ax.get_ylim()[0])/5
    ax.annotate( text, xy=(pow(10,xcenter),ytext), ha='center', va='center', fontsize=fontsize )

ax[2][1].set_ylim(-1,2)
ax[2][1].tick_params(labelleft=False,left=False)
#annotation_line(ax[2][1],1e7,3e8,0,"Radio Waves")
annotation_line(ax[2][1],1e9,3e11,0,"Microwaves")
annotation_line(ax[2][1],3e11,4e14,0,"IR")
annotation_line(ax[2][1],4e14,7.5e14,0,"Visible")
annotation_line(ax[2][1],7.5e14,2e16,0,"UV")
annotation_line(ax[2][1],2e16,1e17,0,"X Rays")
#annotation_line(ax[2][1],2e19,1e23,0,"Gamma Rays")
ax[2,1].set_xlabel(r"$\nu$ ($Hz$)")

ax[2][0].set_ylim(-1,2)
ax[2][0].tick_params(labelleft=False,left=False)
#annotation_line(ax[2][0],1e0,1e3,0,"Radio Waves")
annotation_line(ax[2][0],1e-3,3e-1,0,"Microwaves")
annotation_line(ax[2][0],7.5e-7,1e-3,0,"IR")
annotation_line(ax[2][0],4e-7,7.5e-7,0,"Visible")
annotation_line(ax[2][0],8e-9,4e-7,0,"UV")
annotation_line(ax[2][0],3e-9,8e-9,0,"X Rays")
#annotation_line(ax[2][0],1e-12,1e-11,0,"Gamma Rays")
ax[2,0].set_xlabel(r"$\lambda$ ($m$)")

fig.suptitle(r"Plots for an isotropic synchrotron radiation for $B = 10^{-4}T$, $p=2.0$, $\gamma = 10^{4}$, $\alpha_{p} = 45^{\circ}$.")

plt.savefig("isotropic_sed.png",bbox_inches='tight')

#plt.show()
