#!/usr/bin/python
#-*- coding:utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
from matplotlib.ticker import NullFormatter

A=np.loadtxt("data_boost.res")

nu = A[0,0]
B = A[0,1]
p = A[0,2]

val = []
ydata_I = []

L_ref = A[1,0]
xdata_w = A[2:,0]*180./(np.pi)

for i in range(1,len(A[0,:])):
    val.append(A[1,i])
    ydata_I.append(A[2:,i]/L_ref)

fig,ax = plt.subplots(figsize=(18,9))
plt.subplots_adjust(left=0.06,bottom=0.07,right=0.98,top=0.9,wspace=0.2,hspace=0.3)

for i in range(len(A[0,:])-1):
    ax.plot(xdata_w,ydata_I[i],label=r"$\gamma_{bulk}$"+" = {:0.2f}".format(val[i]))

#ax.set_xscale('log')
ax.set_yscale('log')
ax.yaxis.set_minor_formatter(NullFormatter())
ax.tick_params(top=True, right=True)
ax.set_xlim(0,90)
ax.set_ylim(3e-6,1e7)
#ax.set_ylim(max(ydata_I[-1][ydata_I[-1] != 0])*1e-9,max(ydata_I[-1][ydata_I[-1] != 0])*2)
ax.set_xlabel(r"$\phi$ ($deg$)")
ax.set_ylabel(r"Relative Emissivity")
ax.legend()
ax.set_title("Emissivity for a population of particles")

fig.suptitle(r"Plot for an ordered synchrotron radiation with relativistic Doppler boosting for $\nu=$"+"{:0.2e}Hz, ".format(nu)+r"$p$="+"{}, ".format(p)+r"$B=$"+"{:0.2e}T".format(B))
plt.savefig("boostedflux_phidependency.png",bbox_inches='tight')

#plt.show()
