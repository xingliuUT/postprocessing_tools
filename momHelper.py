#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from momlib import *
from geomHelper import *
from plotHelper import *
from windowFFT import *

zi = complex(0, 1)

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
                   setTime = - 1, \
                   show_plots = False, \
                   plot_format = 'display'):
    momen.set_time(momen.tmom[setTime])
    time = momen.tmom[setTime]
    print 'Reading moments are at t = ', time
    nz = momen.pars['nz0']
    nky = momen.pars['nky0']
    nx = momen.pars['nx0']
    dz = 2.0/nz
    zgrid = np.arange(nz)/float(nz-1)*(2.0-dz)-1.0
    if 'lx_a' in momen.pars:
        xgrid = np.arange(nx)/float(nx-1)*momen.pars['lx_a']+momen.pars['x0']-momen.pars['lx_a']/2.0
    else:
        xgrid = np.arange(nx)/float(nx-1)*momen.pars['lx'] - momen.pars['lx']/2.0
    if zInd == -1 and kyInd != -1 and xInd == -1: 
        deln = momen.dens()[0 : nz, kyInd, 0 : nx]
        tperp = momen.tperp()[0 : nz, kyInd, 0 : nx]
    elif zInd != -1 and kyInd != -1 and  xInd == -1:
        deln = momen.dens()[zInd, kyInd, 0 : nx]
        tperp = momen.tperp()[zInd, kyInd, 0 : nx]
    elif zInd != -1 and kyInd == -1 and  xInd != -1:
        deln = momen.dens()[zInd, 0 : nky, xInd]
        tperp = momen.tperp()[zInd, 0 : nky, xInd]
    deln = deln * momen.pars['rhostar']
    tperp = tperp * momen.pars['rhostar']
    if show_plots:# and i == momen.pars['nky0'] - 1:
        title = 'ky=' + str(kyInd)
        filename = 'n='+str(kyInd*6)+'_dens_tperp_time='+str(np.round(time,4))+'.ps'
#        singlePlot2D(xgrid, zgrid, dens, 'dens_xz', title, filename, 'x', 'z', 'display')
        doublePlot2D(xgrid, zgrid, deln, tperp, 'dens', 'tperp', title, filename, 'x', 'z', plot_format)
    return time, deln, tperp
def momen_xz(momen, \
             geom_coeff, \
             zgrid, \
             kygrid, \
             xgrid, \
             timeInd = -1, \
             show_plots = False, \
             plot_format = 'display'):
    debug = False
    show_raw_plots = False
    q, Cy = q_Cy(geom_coeff)
    nGrid = np.array(kygrid)*momen.pars['n0_global']
    thetaGrid = zgrid * np.pi
    thetaqMatrix = np.outer(thetaGrid, q)
    if debug:
        print 'zi='+str(zi)
        print 'n0='+str(momen.pars['n0_global'])
        print q
        print nGrid
        print thetaGrid
        print thetaqMatrix
    dens_xz = np.zeros((len(zgrid), len(q)),dtype='complex128')
    tperp_xz = np.zeros((len(zgrid), len(q)),dtype='complex128')
    for ky in kygrid:
        time, this_dens, this_tperp = global_moments(momen, -1, ky, -1, timeInd, show_raw_plots, plot_format)
        dens_xz += np.multiply(this_dens, np.exp(zi * nGrid[ky] * thetaqMatrix))
        tperp_xz += np.multiply(this_tperp, np.exp(zi * nGrid[ky] * thetaqMatrix))
        if ky != 0:
            dens_xz += np.multiply(np.conj(this_dens), np.exp(- zi * nGrid[ky] * thetaqMatrix))
            tperp_xz += np.multiply(np.conj(this_tperp), np.exp(- zi * nGrid[ky] * thetaqMatrix))
    if show_plots:
        title = 'time='+str(np.round(time,4))
        filename = 'xz_dens_tperp_time='+str(np.round(time,4))+'.ps'
#        singlePlot2D(xgrid, zgrid, dens_xz, 'dens_xz', title, filename, 'x', 'z', 'display')
        doublePlot2D(xgrid, zgrid, dens_xz, tperp_xz, 'dens_xz', 'tperp_xz', title, filename, 'x', 'z', plot_format)
    return time, dens_xz, tperp_xz

def momen_tx(momen, \
                  geom_coeff, \
                  zgrid, \
                  kygrid, \
                  xgrid, \
                  zInd, \
                  tStart, \
                  tEnd, \
                  show_xz = False, \
                  plot_format = 'display'):
    itStart = np.argmin(abs(np.array(momen.tmom)-tStart))
    itEnd = np.argmin(abs(np.array(momen.tmom)-tEnd))
    tsteps = itEnd - itStart + 1
    tgrid = []
    nz = momen.pars['nz0']
    nx = momen.pars['nx0']
    deln_tx = np.zeros((tsteps, nx),dtype='complex128')
    tperp_tx = np.zeros((tsteps, nx),dtype='complex128')
    for timeInd in range(itStart, itEnd + 1, 350):
        deln_x = np.zeros(nx,dtype='complex128')
        tperp_x = np.zeros(nx,dtype='complex128')
        if show_xz:
            time, dens_xz, tperp_xz = momen_xz(momen, geom_coeff, zgrid, kygrid, xgrid, timeInd, True, plot_format)
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

