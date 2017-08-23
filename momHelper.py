#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from momlib import *
from geomHelper import *
from plotHelper import *
from windowFFT import *

def momen_step_time(momen, \
              show_plots = True):
    momen_time = momen.tmom
    step_time = np.array(momen_time[1:-1]) - np.array(momen_time[0:-2])
    if show_plots:
        plt.plot(step_time)
        plt.title('dens, Tperp')
        plt.ylabel('step time (Lref / cref)')
        plt.show()
def global_moments(momen, \
                   zInd, \
                   kyInd, \
                   xInd, \
                   setTime = - 1):
    momen.set_time(momen.tmom[setTime])
    time = momen.tmom[setTime]
    print 'Reading moments are at t = ', time
    nz = momen.pars['nz0']
    nky = momen.pars['nky0']
    nx = momen.pars['nx0']
    if zInd == -1 and kyInd != -1 and xInd == -1: 
        deln = momen.dens()[0 : nz, kyInd, 0 : nx]
        tperp = momen.tperp()[0 : nz, kyInd, 0 : nx]
    elif zInd != -1 and kyInd != -1 and  xInd == -1:
        deln = momen.dens()[zInd, kyInd, 0 : nx]
        tperp = momen.tperp()[zInd, kyInd, 0 : nx]
    elif zInd != -1 and kyInd == -1 and  xInd != -1:
        deln = momen.dens()[zInd, 0 : nky, xInd]
        tperp = momen.tperp()[zInd, 0 : nky, xInd]
    return time, deln, tperp
def momen_xz(momen, geom_coeff, zgrid, kygrid, xgrid, timeInd = -1, show_plots = False):

    q, Cy = q_Cy(geom_coeff)
    qCy = np.array(q * Cy)
    ymatrix = np.outer(zgrid*np.pi, qCy)
    dens_xz = np.zeros((len(zgrid), len(q)),dtype='complex128')
    tperp_xz = np.zeros((len(zgrid), len(q)),dtype='complex128')
    kygrid = np.array(kygrid) * momen.pars['kymin']
    for i in range(len(kygrid)):
        time, this_dens, this_tperp = global_moments(momen, -1, i, -1, timeInd)
        this_dens = this_dens * momen.pars['rhostar']
        dens_xz += np.multiply(this_dens, np.exp(zi * kygrid[i] * ymatrix))

        this_tperp = this_tperp * momen.pars['rhostar']
        tperp_xz += np.multiply(this_tperp, np.exp(zi * kygrid[i] * ymatrix))
        if i != 0:
            dens_xz += np.multiply(np.conj(this_dens), np.exp(- zi * kygrid[i] * ymatrix))
            tperp_xz += np.multiply(np.conj(this_tperp), np.exp(- zi * kygrid[i] * ymatrix))
        if show_plots:# and i == momen.pars['nky0'] - 1:
            title = 'ky =' + str(i)
            filename = 'tmp.ps'
#            singlePlot2D(xgrid, zgrid, this_dens, 'dens_xz', title, filename, 'x', 'z', 'display')
            doublePlot2D(xgrid, zgrid, this_dens, this_tperp, 'dens_xz', 'tperp_xz', title, filename, 'x', 'z', 'display')
    return time, dens_xz, tperp_xz

def momen_tx(momen, \
                  geom_coeff, \
                  zgrid, \
                  kygrid, \
                  xgrid, \
                  zInd, \
                  tStart, \
                  tEnd, \
                  show_xz = False):
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
        if show_xz:
            time, dens_xz, tperp_xz = momen_xz(momen, geom_coeff, zgrid, kygrid, xgrid, timeInd, True)
        else:
            time, dens_xz, tperp_xz = momen_xz(momen, geom_coeff, zgrid, kygrid, xgrid, timeInd)
        deln_x = dens_xz[zInd,:]
        tperp_x = tperp_xz[zInd,:]
        deln_tx[timeInd - itStart, :] = deln_x.reshape(1, nx)
        tperp_tx[timeInd - itStart, :] = tperp_x.reshape(1, nx)
        tgrid.append(time)
    return tgrid, deln_tx, tperp_tx

def momen_fx(field_tx, tgrid, nf, lf):
    tsteps, nx = np.shape(field_tx)
    field_fx = np.zeros((nf, nx), dtype='complex128')
    for i in range(nx):
        fgrid, field_f = windowFFT(np.array(tgrid), field_tx[:,i], nf, lf, str(i/float(nx)))
        field_f = field_f.reshape(nf)
        field_fx[:,i] = field_f
    return fgrid, field_fx
def momen_tky(momen, \
                  zInd, \
                  kyInd, \
                  xInd, \
                  tStart, \
                  tEnd):
    nky = momen.pars['nky0']
    itStart = np.argmin(abs(np.array(momen.tmom)-tStart))
    itEnd = np.argmin(abs(np.array(momen.tmom)-tEnd))
    tsteps = itEnd - itStart + 1
    tgrid = []
    deln_tky = np.zeros((tsteps, nky),dtype='complex128')
    tperp_tky = np.zeros((tsteps, nky),dtype='complex128')
    for timeInd in range(itStart, itEnd + 1):
        time, this_deln, this_tperp = global_moments(momen, zInd, kyInd, xInd, timeInd)
        deln_tky[timeInd - itStart, :] = this_deln.reshape(1, nky)
        tperp_tky[timeInd - itStart, :] = this_tperp.reshape(1, nky)
        tgrid.append(time)
    return tgrid, deln_tky, tperp_tky

