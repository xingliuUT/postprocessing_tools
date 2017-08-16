#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from fieldlib import *

def global_eigenfunctions(field, \
                          zInd = -1, \
                          kyInd = 0, \
                          xInd = -1, \
                          setTime = -1):
    field.set_time(field.tfld[setTime])
    time = field.tfld[setTime]
    print 'Reading eigenfunctions are at t = ', time
    nz = field.nz
    nx = field.nx
    phi = np.zeros((nz, nx), dtype='complex128')
    apar = np.zeros((nz, nx), dtype='complex128')
    if zInd == -1 and xInd == -1:
        phi = field.phi()[0 : nz, kyInd, 0 : nx]
        apar = field.apar()[0 : nz, kyInd, 0 : nx]
    elif xInd == -1:
        phi = field.phi()[zInd, kyInd, 0 : nx]
        apar = field.apar()[zInd, kyInd, 0 : nx]
    return time, phi, apar

def t_avg_global_eigenfns(field, \
                          kyInd = 0, \
                          tStart = 0, \
                          tEnd = 0):
    itStart = np.argmin(abs(np.array(field.tfld) - tStart))
    itEnd = np.argmin(abs(np.array(field.tfld) - tEnd))
    phi = np.zeros((field.nz, field.nx), dtype='complex128')
    apar = np.zeros((field.nz, field.nx), dtype='complex128')
    for timeInd in range(itStart, itEnd + 1):
        t, this_phi, this_apar = global_eigenfunctions(field, \
                                 kyInd, timeInd)
        phi += this_phi
        apar += this_apar
    return phi / (itEnd - itStart + 1), apar / (itEnd - itStart + 1)

def field_step_time(field, \
              show_plots = True):
    field_time = field.tfld
    step_time = np.array(field_time[1:-1]) - np.array(field_time[0:-2])
    if show_plots:
        plt.plot(step_time)
        plt.title('Phi and Apar')
        plt.ylabel('step time (Lref / cref)')
        plt.show()

def reflectometer(field, \
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
    print len(tgrid)
    print np.shape(phi_tx)
    print np.shape(apar_tx)
    return tgrid, phi_tx, apar_tx
