#!/usr/bin/python
#-*- coding:utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d

def annotation_line(ax, xmin, xmax, y, text, ytext=0, linecolor='black', linewidth=1, fontsize=8):
    ax.annotate('', xy=(xmin, y), xytext=(xmax, y), xycoords='data', textcoords='data',arrowprops={'arrowstyle': '|-|', 'color':linecolor, 'linewidth':linewidth})
    xcenter = np.log10(xmin)+(np.log10(xmax)-np.log10(xmin))/2
    if ytext==0:
        ytext = y+(ax.get_ylim()[1]-ax.get_ylim()[0])/5
    ax.annotate( text, xy=(pow(10,xcenter),ytext), ha='center', va='center', fontsize=fontsize )

A=np.loadtxt("data_ord_yp.res")

val = []
ydata_I = []

y_p = A[0,0]
xdata_w = A[1:,0]
xdata_l = 3e8/xdata_w

for i in range(1,len(A[0,:])):
    val.append(A[0,i])
    ydata_I.append(A[1:,i])

fig,ax = plt.subplots(2,2,sharex='col',gridspec_kw={'height_ratios': [3,0.20]},figsize=(18,9))
plt.subplots_adjust(left=0.06,bottom=0.07,right=0.98,top=0.9,wspace=0.2,hspace=0)

for i in range(len(A[0,:])-1):
    #ax[0][0].plot(xdata_w,ydata_I[i],label=r"$\phi$"+" = {:0.2f}deg".format(val[i]*180./np.pi))
    #ax[0][0].plot(xdata_w,ydata_I[i],label=r"$B$"+" = {:0.2e}T".format(val[i]))
    ax[0][0].plot(xdata_w,ydata_I[i],label=r"$\gamma_{p}$"+" = {:0.2e}".format(val[i]))
    #ax[0][0].plot(xdata_w,ydata_I[i],label=r"$p$"+" = {:0.2f}".format(val[i]))

ax[0][0].set_xscale('log')
ax[0][0].set_yscale('log')
ax[0][0].set_xlim(1e7,1e23)
#ax[0][0].set_ylim(1e-42,1e-28)
ax[0][0].set_ylim(max(ydata_I[0][ydata_I[0] != 0])*1e-28,max(ydata_I[0][ydata_I[0] != 0])*10)
ax[0][0].set_xlabel(r"$\nu$ ($Hz$)")
ax[0][0].set_ylabel(r"Emissivity ($W \cdot Hz^{-1}$)")
ax[0][0].legend()
ax[0][0].set_title("Emissivity for a population of particles")

for i in range(len(A[0,:])-1):
    #ax[0][1].plot(xdata_l,ydata_I[i]*3e8/xdata_l**2,label=r"$\phi$"+" = {:0.2f}deg".format(val[i]*180./np.pi))
    #ax[0][1].plot(xdata_l,ydata_I[i]*3e8/xdata_l**2,label=r"$B$"+" = {:0.2e}T".format(val[i]))
    ax[0][1].plot(xdata_l,ydata_I[i]*3e8/xdata_l**2,label=r"$\gamma_{p}$"+" = {:0.2e}".format(val[i]))
    #ax[0][1].plot(xdata_l,ydata_I[i]*3e8/xdata_l**2,label="p = {:0.2f}".format(val[i]))

ax[0][1].set_xscale('log')
ax[0][1].set_yscale('log')
ax[0][1].set_xlabel(r"$\lambda$ ($m$)")
ax[0][1].set_ylabel(r"Emissivity ($W \cdot m^{-1}$)")
ax[0][1].set_xlim(1e-12,1e3)
#ax[0][1].set_ylim(1e-27,1e-9)
ax[0][1].set_ylim(max(ydata_I[-1]*3e8/xdata_l**2)*1e-28,max(ydata_I[-1]*3e8/xdata_l**2)*10)
ax[0][1].legend()
ax[0][1].set_title("Emissivity for a population of particles")

ax[1][0].set_ylim(-1,2)
ax[1][0].tick_params(labelleft=False,left=False)
annotation_line(ax[1][0],1e7,3e8,0,"Radio Waves")
annotation_line(ax[1][0],3e8,3e11,0,"Microwaves")
annotation_line(ax[1][0],3e11,4e14,0,"IR")
annotation_line(ax[1][0],4e14,7.5e14,0,"Visible")
annotation_line(ax[1][0],7.5e14,2e16,0,"UV")
annotation_line(ax[1][0],2e16,2e19,0,"X Rays")
annotation_line(ax[1][0],2e19,1e23,0,"Gamma Rays")
ax[1][0].set_xlabel(r"$\nu$ ($Hz$)")

ax[1][1].set_ylim(-1,2)
ax[1][1].tick_params(labelleft=False,left=False)
annotation_line(ax[1][1],1e0,1e3,0,"Radio Waves")
annotation_line(ax[1][1],1e-3,1e0,0,"Microwaves")
annotation_line(ax[1][1],7.5e-7,1e-3,0,"IR")
annotation_line(ax[1][1],4e-7,7.5e-7,0,"Visible")
annotation_line(ax[1][1],8e-9,4e-7,0,"UV")
annotation_line(ax[1][1],1e-11,8e-9,0,"X Rays")
annotation_line(ax[1][1],1e-12,1e-11,0,"Gamma Rays")
ax[1][1].set_xlabel(r"$\lambda$ ($m$)")

fig.suptitle(r"Plots for an ordered synchrotron radiation.\n with $\phi=45$, $\gamma_{p}=10^{4.5}$, $B=10^{-4}T$, $p=2.0$ (when not changing)")
plt.savefig("ordered_sed.png",bbox_inches='tight')

#plt.show()
