#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from fieldlib import *

def global_eigenfunctions(pars,suffix,show_plots,setTime=-1):

    field = fieldfile('field'+suffix,pars)

    if (setTime == -1):
        field.set_time(field.tfld[setTime])
        print 'Reading eigenfunctions are at t = ', field.tfld[setTime]
    else:
        isetTime = np.argmin(abs(np.array(field.tfld)-setTime))
        field.set_time(field.tfld[isetTime])
        print 'Reading eigenfunctions are at t = ', field.tfld[isetTime]

    phi = np.zeros((field.nz, field.nx), dtype='complex128')
    apar = np.zeros((field.nz, field.nx), dtype='complex128')

    phi = field.phi()[:,0,:]
    apar = field.apar()[:,0,:]

    dz = 2.0/field.nz
    zgrid = np.arange(field.nz)/float(field.nz-1)*(2.0-dz)-1.0
    if 'lx_a' in pars:
        xgrid = np.arange(field.nx)/float(field.nx-1)*pars['lx_a']+pars['x0']-pars['lx_a']/2.0
    else:
        xgrid = np.arange(field.nx)/float(field.nx-1)*pars['lx'] - pars['lx']/2.0

    if show_plots:
        plt.figure(figsize=(6.,8.))
        fig=plt.gcf()
        plt.subplot(3,1,1)
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.title(r'$|\phi|$')
        plt.contourf(xgrid,zgrid,np.abs(phi),70)
        cb1=plt.colorbar()
        plt.subplot(3,1,2)
        plt.title(r'$Re[\phi]$')
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.contourf(xgrid,zgrid,np.real(phi),70)
        plt.colorbar()
        plt.subplot(3,1,3)
        plt.title(r'$Im[\phi]$')
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.contourf(xgrid,zgrid,np.imag(phi),70)
        plt.colorbar()
        plt.tight_layout()
        plt.show()
    if show_plots:
        plt.figure(figsize=(6.,8.))
        fig=plt.gcf()
        plt.subplot(3,1,1)
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.title(r'$|A_{||}|$')
        plt.contourf(xgrid,zgrid,np.abs(apar),70)
        cb1=plt.colorbar()
        plt.subplot(3,1,2)
        plt.title(r'$Re[A_{||}]$')
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.contourf(xgrid,zgrid,np.real(apar),70)
        plt.colorbar()
        plt.subplot(3,1,3)
        plt.title(r'$Im[A_{||}]$')
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.contourf(xgrid,zgrid,np.imag(apar),70)
        plt.colorbar()
        plt.tight_layout()
        plt.show()

    return phi, apar
