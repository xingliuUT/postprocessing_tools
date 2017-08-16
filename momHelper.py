#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from momlib import *

def global_moments(momen, \
                   zInd, \
                   kyInd, \
                   xInd, \
                   setTime = - 1):
    momen.set_time(momen.tmom[setTime])
    time = momen.tmom[setTime]
    print 'Reading moments are at t = ', time
    nz = momen.pars['nz0']
    nx = momen.pars['nx0']
    deln = np.zeros((nz,nx),dtype='complex128')
    tperp = np.zeros((nz,nx),dtype='complex128')
    if zInd == -1 and xInd == -1: 
        deln = momen.dens()[0 : nz, kyInd, 0 : nx]
        tperp = momen.tperp()[0 : nz, kyInd, 0 : nx]
    elif xInd == -1:
        deln = momen.dens()[zInd, kyInd, 0 : nx]
        tperp = momen.tperp()[zInd, kyInd, 0 : nx]
    return time, deln, tperp

def t_avg_global_moms(momen, \
                      kyInd = 0, \
                      tStart = 0., \
                      tEnd = 0.):

    itStart = np.argmin(abs(np.array(momen.tmom)-tStart))
    itEnd = np.argmin(abs(np.array(momen.tmom)-tEnd))

    nz = momen.pars['nz0']
    nx = momen.pars['nx0']

    deln = np.zeros((nz,nx),dtype='complex128')
    tperp = np.zeros((nz,nx),dtype='complex128')

    for timeInd in range(itStart, itEnd + 1):
        t, this_deln, this_tperp = global_moments(momen, kyInd, timeInd)
        deln += this_deln
        tperp += this_tperp

    return deln / (itEnd - itStart + 1), tperp / (itEnd - itStart + 1)

def momen_step_time(momen, \
              show_plots = True):
    momen_time = momen.tmom
    step_time = np.array(momen_time[1:-1]) - np.array(momen_time[0:-2])
    if show_plots:
        plt.plot(step_time)
        plt.title('dens, Tperp')
        plt.ylabel('step time (Lref / cref)')
        plt.show()
def momen_reflectometer(momen, \
                  zInd, \
                  kygrid, \
                  xInd, \
                  tStart, \
                  tEnd):

    itStart = np.argmin(abs(np.array(momen.tmom)-tStart))
    itEnd = np.argmin(abs(np.array(momen.tmom)-tEnd))
    tsteps = itEnd - itStart + 1
    tgrid = []
    nz = momen.pars['nz0']
    nx = momen.pars['nx0']
    deln_tx = np.zeros((tsteps, nx),dtype='complex128')
    tperp_tx = np.zeros((tsteps, nx),dtype='complex128')

    for timeInd in range(itStart, itEnd + 1):
        deln_x = np.zeros(nx,dtype='complex128')
        tperp_x = np.zeros(nx,dtype='complex128')
        for ky in kygrid:
            time, this_deln, this_tperp = global_moments(momen, zInd, ky, xInd, timeInd)
            deln_x += this_deln
            tperp_x += this_tperp
        deln_tx[timeInd - itStart, :] = deln_x.reshape(1, nx)
        tperp_tx[timeInd - itStart, :] = tperp_x.reshape(1, nx)
        tgrid.append(time)
    return tgrid, deln_tx, tperp_tx
