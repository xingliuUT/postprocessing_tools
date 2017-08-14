#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from fieldlib import *

def global_eigenfunctions(pars, \
                          suffix, \
                          kyInd = 0, \
                          show_plots = True, \
                          setTime=-1):

    field = fieldfile('field'+suffix,pars)

    if (setTime == -1):
        field.set_time(field.tfld[setTime])
        time = field.tfld[setTime]
    else:
        isetTime = np.argmin(abs(np.array(field.tfld)-setTime))
        field.set_time(field.tfld[isetTime])
        time = field.tfld[isetTime]
    print 'Reading eigenfunctions are at t = ', time

    phi = np.zeros((field.nz, field.nx), dtype='complex128')
    apar = np.zeros((field.nz, field.nx), dtype='complex128')

    phi = field.phi()[:, kyInd,:]
    apar = field.apar()[:, kyInd,:]

    dz = 2.0/field.nz
    zgrid = np.arange(field.nz)/float(field.nz-1)*(2.0-dz)-1.0
    if 'lx_a' in pars:
        xgrid = np.arange(field.nx)/float(field.nx-1)*pars['lx_a']+pars['x0']-pars['lx_a']/2.0
    else:
        xgrid = np.arange(field.nx)/float(field.nx-1)*pars['lx'] - pars['lx']/2.0

    if show_plots:
        plt.figure(figsize=(12.,8.))
        fig=plt.gcf()
        plt.subplot(3,2,1)
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.title(r'$|\phi|$')
        plt.contourf(xgrid,zgrid,np.abs(phi),70)
        plt.colorbar()
        plt.subplot(3,2,3)
        plt.title(r'$Re[\phi]$')
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.contourf(xgrid,zgrid,np.real(phi),70)
        plt.colorbar()
        plt.subplot(3,2,5)
        plt.title(r'$Im[\phi]$')
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.contourf(xgrid,zgrid,np.imag(phi),70)
        plt.colorbar()
        plt.subplot(3,2,2)
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.title(r'$|A_{||}|$')
        plt.contourf(xgrid,zgrid,np.abs(apar),70)
        plt.colorbar()
        plt.subplot(3,2,4)
        plt.title(r'$Re[A_{||}]$')
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.contourf(xgrid,zgrid,np.real(apar),70)
        plt.colorbar()
        plt.subplot(3,2,6)
        plt.title(r'$Im[A_{||}]$')
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.contourf(xgrid,zgrid,np.imag(apar),70)
        plt.colorbar()
        plt.tight_layout()
        plt.suptitle('time =' + str(np.round(time,1)) + ', ky index =' + str(kyInd))
        plt.show()

    return phi, apar