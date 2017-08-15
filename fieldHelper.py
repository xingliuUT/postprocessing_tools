#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from fieldlib import *


def global_eigenfunctions(field, \
                          kyInd = 0, \
                          setTime = -1):

#    if (setTime == -1):
    field.set_time(field.tfld[setTime])
    time = field.tfld[setTime]
#    else:
#        isetTime = np.argmin(abs(np.array(field.tfld)-setTime))
#        field.set_time(field.tfld[isetTime])
#        time = field.tfld[isetTime]
    print 'Reading eigenfunctions are at t = ', time

    phi = np.zeros((field.nz, field.nx), dtype='complex128')
    apar = np.zeros((field.nz, field.nx), dtype='complex128')

    phi = field.phi()[:, kyInd,:]
    apar = field.apar()[:, kyInd,:]

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
