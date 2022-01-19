#!/usr/bin/python
#-*- coding:utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

def annotation_line(ax, xmin, xmax, y, text, ytext=0, linecolor='black', linewidth=1, fontsize=8):
    ax.annotate('', xy=(xmin, y), xytext=(xmax, y), xycoords='data', textcoords='data',arrowprops={'arrowstyle': '|-|', 'color':linecolor, 'linewidth':linewidth})
    #ax.annotate('', xy=(xmin, y), xytext=(xmax, y), xycoords='data', textcoords='data',arrowprops={'arrowstyle': '<->', 'color':linecolor, 'linewidth':linewidth})
    xcenter = np.log10(xmin)+(np.log10(xmax)-np.log10(xmin))/2
    if ytext==0:
        ytext = y+(ax.get_ylim()[1]-ax.get_ylim()[0])/5
    ax.annotate( text, xy=(pow(10,xcenter),ytext), ha='center', va='center', fontsize=fontsize )

A = np.loadtxt('data_pola.res')

gam = A[0,0]
nu = A[1:,0]

val_test = []

para = []
orth = []
d = []

for i in range(int(len(A[0,:])/3)):
    val_test.append(A[0,1+3*i]*180./np.pi)
    para.append(A[1:,1+3*i])
    orth.append(A[1:,2+3*i])
    d.append(A[1:,3+3*i])

fig,ax = plt.subplots(3,1,sharex=True,gridspec_kw={'height_ratios': [3,1,0.25]},figsize=(18,9))
plt.subplots_adjust(left=0.10,bottom=0.10,right=0.90,top=0.90,wspace=0.5,hspace=0)

for i in range(len(val_test)):
    ax[0].plot(nu,para[i],label=r"$I_{\parallel}$, $\phi=$"+"{:0.2f}deg".format(val_test[i]))
    ax[0].plot(nu,orth[i],label=r"$I_{\perp}$, $\phi=$"+"{:0.2f}deg".format(val_test[i]))
    ax[1].plot(nu,d[i]*100.,label=r"$\phi=$"+"{:0.2f}deg".format(val_test[i]))

ax[0].set_xscale('log')
ax[0].set_xlim(1e7,1e15)
ax[0].set_ylabel(r"Intensity of polarisation components ($W \cdot Hz^{-1}$)")
ax[0].legend()

ax[1].set_xscale('log')
ax[1].set_xlim(1e7,1e15)
ax[1].set_ylim(0,100.)
ax[1].set_ylabel('Polarisation power')
ax[1].set_xlabel('Frequency (Hz)')
ax[1].legend()

ax[2].set_ylim(-1,2)
ax[2].tick_params(labelleft=False,left=False)
annotation_line(ax[2],1e7,3e8,0,"Radio Waves")
annotation_line(ax[2],3e8,3e11,0,"Microwaves")
annotation_line(ax[2],3e11,4e14,0,"IR")
annotation_line(ax[2],4e14,7.5e14,0,"Visible")
ax[2].set_xlabel('Frequency (Hz)')

fig.suptitle(r"Synchrotron emission polarisation as a function of frequency for different values of observation angle $\phi$")
plt.savefig("synchrotron_polarisation.png",bbox_inches='tight')

#plt.show()


