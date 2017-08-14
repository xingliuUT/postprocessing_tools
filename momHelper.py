#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from momlib import *

def global_moments(pars,suffix,kyInd,show_plots,setTime=-1):

    momen = momfile('mom_e'+suffix,pars)
    if (setTime == -1):
        momen.set_time(momen.tmom[setTime])
        time = momen.tmom[setTime]
    else:
        isetTime = np.argmin(abs(np.array(momen.tmom)-setTime))
        momen.set_time(momen.tmom[isetTime])
        time = momen.tmom[isetTime]
    print 'Reading momentss are at t = ', time

    nz = pars['nz0']
    nx = pars['nx0']

    upar = np.zeros((nz,nx),dtype='complex128')
    deln = np.zeros((nz,nx),dtype='complex128')
    tpar = np.zeros((nz,nx),dtype='complex128')
    tperp = np.zeros((nz,nx),dtype='complex128')
    qpar = np.zeros((nz,nx),dtype='complex128')
    qperp = np.zeros((nz,nx),dtype='complex128')

    
    deln = momen.dens()[:,kyInd,:]
    tperp = momen.tperp()[:,kyInd,:]

    dz = 2.0/nz
    zgrid = np.arange(nz)/float(nz-1)*(2.0-dz)-1.0
    if 'lx_a' in pars:
        xgrid = np.arange(nx)/float(nx-1)*pars['lx_a']+pars['x0']-pars['lx_a']/2.0
    else:
        xgrid = np.arange(nx)/float(nx-1)*pars['lx'] - pars['lx']/2.0

    if show_plots:
        plt.figure(figsize=(12.,8.))
        fig=plt.gcf()
        plt.subplot(3,2,1)
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.title(r'$|T_{\perp}|$')
        plt.contourf(xgrid,zgrid,np.abs(tperp),70)
        plt.colorbar()
        plt.subplot(3,2,3)
        plt.title(r'$Re[T_{\perp}]$')
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.contourf(xgrid,zgrid,np.real(tperp),70)
        plt.colorbar()
        plt.subplot(3,2,5)
        plt.title(r'$Im[T_{\perp}]$')
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.contourf(xgrid,zgrid,np.imag(tperp),70)
        plt.colorbar()
        plt.subplot(3,2,2)
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.title(r'$|n|$')
        plt.contourf(xgrid,zgrid,np.abs(deln),70)
        plt.colorbar()
        plt.subplot(3,2,4)
        plt.title(r'$Re[n]$')
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.contourf(xgrid,zgrid,np.real(deln),70)
        plt.colorbar()
        plt.subplot(3,2,6)
        plt.title(r'$Im[n]$')
        plt.ylabel(r'$z/\pi$',fontsize=13)
        plt.xlabel(r'$\rho_{tor}$',fontsize=13)
        plt.contourf(xgrid,zgrid,np.imag(deln),70)
        plt.colorbar()
        plt.tight_layout()
        plt.suptitle('time =' + str(np.round(time,1)) + ', ky index =' + str(kyInd))
        plt.show()

    return deln,tperp

