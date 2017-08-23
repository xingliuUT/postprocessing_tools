#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from fieldlib import *
from geomHelper import *
from plotHelper import *

zi = complex(0,1)
def field_step_time(field, \
                    show_plots = True):
    field_time = field.tfld
    step_time = np.array(field_time[1:-1]) - np.array(field_time[0:-2])
    if show_plots:
        plt.plot(step_time)
        plt.title('Phi and Apar')
        plt.ylabel('step time (Lref / cref)')
        plt.show()
def global_eigenfunctions(field, \
                          zInd, \
                          kyInd, \
                          xInd, \
                          setTime = -1, \
                          show_plots = False, \
                          plot_format = 'display'):
    field.set_time(field.tfld[setTime])
    time = field.tfld[setTime]
    print 'Reading eigenfunctions are at t = ', time
    nz = field.nz
    ny = field.ny
    nx = field.nx
    dz = 2.0/nz
    zgrid = np.arange(nz)/float(nz-1)*(2.0-dz)-1.0
    if 'lx_a' in field.pars:
        xgrid = np.arange(nx)/float(nx-1)*field.pars['lx_a']+field.pars['x0']-field.pars['lx_a']/2.0
    else:
        xgrid = np.arange(nx)/float(nx-1)*field.pars['lx'] - field.pars['lx']/2.0
    if zInd == -1 and kyInd != -1 and xInd == -1:
        phi = field.phi()[0 : nz, kyInd, 0 : nx]
        apar = field.apar()[0 : nz, kyInd, 0 : nx]
    elif zInd != -1 and kyInd != -1 and xInd == -1:
        phi = field.phi()[zInd, kyInd, 0 : nx]
        apar = field.apar()[zInd, kyInd, 0 : nx]
    elif zInd != -1 and kyInd == -1 and xInd != -1:
        phi = field.phi()[zInd, 0:ny, xInd]
        apar = field.apar()[zInd, 0:ny, xInd]
    phi = phi*field.pars['rhostar']
    apar = apar*field.pars['rhostar']
    if show_plots:
        title = 'ky=' + str(kyInd)
        filename = 'n='+str(kyInd*6)+'_phi_apar_time='+str(np.round(time,4))+'.ps'
        doublePlot2D(xgrid, zgrid, phi, apar, 'phi', 'apar', title, filename, 'x', 'z', plot_format)
    return time, phi, apar
def field_xz(field, \
             geom_coeff, \
             zgrid, \
             kygrid, \
             xgrid, \
             timeInd = -1, \
             show_plots = False, \
             plot_format = 'display'):
    show_raw_plots = False
    q, Cy = q_Cy(geom_coeff)
    nGrid = np.array(kygrid)*field.pars['n0_global']
    thetaGrid = zgrid * np.pi
    thetaqMatrix = np.outer(thetaGrid, q)
    phi_xz = np.zeros((len(zgrid),len(q)), dtype = 'complex128')
    apar_xz = np.zeros((len(zgrid),len(q)), dtype = 'complex128')
    for ky in kygrid:
        time, this_phi, this_apar = global_eigenfunctions(field, -1, ky, -1, timeInd, show_raw_plots, plot_format)
        phi_xz += np.multiply(this_phi, np.exp(zi * nGrid[ky] * thetaqMatrix))
        apar_xz += np.multiply(this_apar, np.exp(zi * nGrid[ky] * thetaqMatrix))
        if ky != 0:
            phi_xz += np.multiply(np.conj(this_phi), np.exp(- zi * nGrid[ky] * thetaqMatrix))
            apar_xz += np.multiply(np.conj(this_apar), np.exp(- zi * nGrid[ky] * thetaqMatrix))
    if show_plots:
        title = 'time='+str(np.round(time,4))
        filename = 'phi_apar_time='+str(np.round(time,4))+'.ps'
#        singlePlot2D(xgrid, zgrid, phi_xz, 'dens_xz', title, filename, 'x', 'z', 'display')
        doublePlot2D(xgrid, zgrid, phi_xz, apar_xz, 'phi_xz', 'apar_xz', title, filename, 'x', 'z', plot_format)
    return time, phi_xz, apar_xz
            

def field_tx(field, \
                  zInd, \
                  kygrid, \
                  xInd, \
                  tStart, \
                  tEnd):

    itStart = np.argmin(abs(np.array(field.tfld) - tStart))
    itEnd = np.argmin(abs(np.array(field.tfld) - tEnd))
    tsteps = itEnd - itStart + 1
    tgrid = []
    phi_tx = np.zeros((tsteps, field.nx), dtype='complex128')
    apar_tx = np.zeros((tsteps, field.nx), dtype='complex128')
    for timeInd in range(itStart, itEnd + 1):
        phi_x = np.zeros(field.nx, dtype='complex128')
        apar_x = np.zeros(field.nx, dtype='complex128')
        for ky in kygrid:
            time, this_phi, this_apar = global_eigenfunctions(field, zInd, ky, xInd, timeInd)
            phi_x += this_phi
            apar_x += this_apar
        phi_tx[timeInd - itStart, :] = phi_x.reshape(1, field.nx)
        apar_tx[timeInd - itStart, :] = apar_x.reshape(1, field.nx)
        tgrid.append(time)
    #print len(tgrid)
    #print np.shape(phi_tx)
    #print np.shape(apar_tx)
    return tgrid, phi_tx, apar_tx
