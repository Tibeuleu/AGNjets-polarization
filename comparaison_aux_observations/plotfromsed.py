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

z = 0.77
Dl = 3407.27e6*3.0857e16

SED = np.loadtxt('photandseds.tbl',comments='|')
INT = np.loadtxt('Flux_modele.res')
MOD = np.loadtxt('data_ord_p_3.res')

n_val = len(MOD[0,:])
n_point = len(SED[:,0])
n_ratio = int(float(n_point)/2.5)

X_sed = SED[:,0]
Y_sed = SED[:,1]*1e-26
yerr_sed = SED[:,2:3]*1e-26

X_mod = MOD[1:,0]
Y_mod = []
val_test = []

X_int = INT[1:,0]
Y_int = []

for i in range(1,n_val):

    val_test.append(MOD[0,i])

    ratio = Y_sed[n_ratio]/(INT[n_ratio+1,i]*pow(1+z,(val_test[i-1]+1.)/2.)/(4.*np.pi*Dl*Dl))
    Y_mod.append(ratio*MOD[1:,i]*pow(1+z,(val_test[i-1]+1.)/2.)/(4.*np.pi*Dl*Dl))

    Y_int.append(ratio*INT[1:,i]*pow(1+z,(val_test[i-1]+1.)/2.)/(4.*np.pi*Dl*Dl))

fig,ax = plt.subplots(3,1,sharex=True,gridspec_kw={'height_ratios': [3,1,0.25]},figsize=(18,9))
plt.subplots_adjust(left=0.10,bottom=0.10,right=0.90,top=0.85,wspace=0.5,hspace=0)

ax[0].errorbar(X_sed,Y_sed,yerr_sed,fmt="ok")
for i in range(n_val-1):
    ax[0].plot(X_mod,Y_mod[i],'--',color="C{}".format(i),label="p={}".format(val_test[i]))
    ax[0].plot(X_int,Y_int[i],'+',color="C{}".format(i))
ax[0].set_yscale('log')
ax[0].set_xscale('log')
ax[0].set_xlim(1e7,1e15)
ax[0].set_ylim(min(Y_sed[Y_sed != 0])*1e-3,max(Y_sed[Y_sed != 0])*1e1)
ax[0].set_ylabel(r"Flux density ($W \cdot m^{-2} \cdot Hz^{-1}$)")
ax[0].legend()

for i in range(n_val-1):
    ax[1].errorbar(X_sed,(Y_sed/Y_int[i]),yerr_sed,fmt='o',color="C{}".format(i),label="p={}".format(val_test[i]))
ax[1].plot([1e7,1e15],[1,1],'--k')
ax[1].set_xscale('log')
ax[1].set_yscale('log')
ax[1].set_xlim(1e7,1e15)
ax[1].set_ylabel('Ratio')
ax[1].set_xlabel('Frequency (Hz)')

ax[2].set_ylim(-1,2)
ax[2].tick_params(labelleft=False,left=False)
annotation_line(ax[2],1e7,3e8,0,"Radio Waves")
annotation_line(ax[2],3e8,3e11,0,"Microwaves")
annotation_line(ax[2],3e11,4e14,0,"IR")
annotation_line(ax[2],4e14,7.5e14,0,"Visible")
ax[2].set_xlabel('Frequency (Hz)')

fig.suptitle(r"3C175 spectral energy distribution fit using a toy model and calibrating on photon index $p$."+"\n (SED source : NED/NASA)")

plt.savefig("SED_fit.png",bbox_inches='tight')

#plt.show()
