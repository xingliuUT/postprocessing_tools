#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from momlib import *

def global_moments(momen, \
                   kyInd, \
                   setTime = - 1):

    #if (setTime == -1):
    momen.set_time(momen.tmom[setTime])
    time = momen.tmom[setTime]
    #else:
    #    isetTime = np.argmin(abs(np.array(momen.tmom)-setTime))
    #    momen.set_time(momen.tmom[isetTime])
    #    time = momen.tmom[isetTime]
    print 'Reading momentss are at t = ', time

    nz = momen.pars['nz0']
    nx = momen.pars['nx0']

    deln = np.zeros((nz,nx),dtype='complex128')
    tperp = np.zeros((nz,nx),dtype='complex128')

    
    deln = momen.dens()[:,kyInd,:]
    tperp = momen.tperp()[:,kyInd,:]

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

