#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from fieldHelper import *
from momHelper import *
from parIOHelper import *
from geomHelper import *
from plotHelper import *
from windowFFT import *
import sys

suffix = sys.argv[1]

if not suffix =='.dat':
   suffix = '_'+suffix

tStart = float(sys.argv[2])
tEnd = float(sys.argv[3])

pars = init_read_parameters(suffix)
momen = momfile('mom_e'+suffix,pars)

zi = complex(0, 1)
nz = pars['nz0']
nx = pars['nx0']
ny = pars['nky0']
dz = 2.0/nz
zgrid = np.arange(nz)/float(nz-1)*(2.0-dz)-1.0
if 'lx_a' in pars:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx_a']+pars['x0']-pars['lx_a']/2.0
else:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx'] - pars['lx']/2.0

show_plots = True
plot_format = 'display'
#plot_format = 'ps'

kygrid = range(pars['nky0'])
#kygrid = [0]
zInd = nz / 2
xInd = nx * 5 / 8

nf = 200
lf = 10.

#for kyInd in kygrid:
if 1 == 1:
    kyInd = -1
    tgrid, dens_tky, tperp_tky = momen_tky(momen, \
                  zInd, \
                  kyInd, \
                  xInd, \
                  tStart, \
                  tEnd)
    title = ' '
    filename = 'dens_tky01.ps'

    doublePlot2D(kygrid, tgrid, dens_tky, tperp_tky, 'dens_tky', 'tperp_tky', title, filename, 'ky', 't',plot_format)
if 1 == 1:
#    plot_format = 'ps'
    dens_fky = np.empty((nf,ny),dtype = 'complex128')
    tperp_fky = np.empty((nf,ny),dtype = 'complex128')
    show_plots = False
    for ky in kygrid:
        fgrid, dens_f = windowFFT(tgrid, dens_tky[:,ky], nf, lf, 'dens_ky=' + str(ky), show_plots, plot_format)
        fgrid, tperp_f = windowFFT(tgrid, tperp_tky[:,ky], nf, lf, 'tperp_ky=' + str(ky), show_plots, plot_format)
        dens_f = np.minimum(dens_f, 0.002)
        tperp_f = np.minimum(tperp_f, 0.002)
        dens_fky[:,ky] = dens_f.reshape(1, nf)
        tperp_fky[:,ky] = tperp_f.reshape(1, nf)
    doublePlot2D(kygrid, fgrid, dens_fky, tperp_fky, 'dens_fky', 'tperp_fky', title, filename, 'ky', 'f',plot_format)

