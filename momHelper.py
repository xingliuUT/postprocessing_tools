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

    return time, deln,tperp

