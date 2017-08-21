#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from fieldHelper import *
from momHelper import *
from parIOHelper import *
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
#momen_step_time(momen)

nz = pars['nz0']
nx = pars['nx0']
if 'lx_a' in pars:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx_a']+pars['x0']-pars['lx_a']/2.0
else:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx'] - pars['lx']/2.0

#plot_format = 'display'
plot_format = 'ps'
nf = 400
lf = 20.
minVal = 1.

#kygrid = np.array(range(pars['nky0'])) * pars['kymin']
kygrid = range(pars['nky0'])

zInd = nz/2
kyInd = -1
#xInd = nx / 2
for xInd in range(nx/2, nx * 3 / 4, 8):
    print xgrid[xInd]
    tgrid, dens_tky, tperp_tky = momen_tky(momen,zInd,kyInd,xInd,tStart,tEnd)
    if 1 == 0:
        title = 'xgrid: time, ygrid: ky, x ='+str(xgrid[xInd])
        filename = 'tky_dens_tperp.ps'
        doublePlot2D(kygrid, tgrid, dens_tky, tperp_tky, 'dens', 'tperp', title, filename, 'ky', 't', plot_format)

    fgrid, dens_fky = momen_fx(dens_tky, tgrid, nf, lf)
#    fgrid, tperp_fky = momen_fx(tperp_tky, tgrid, nf, lf)
    if 1 == 1:
        title = ' '
        filename = 'fky_dens_x='+str(xgrid[xInd])+'.ps'
#        doublePlot2D(kygrid, fgrid, dens_fky, tperp_fky, 'dens', 'tperp', title, filename, 'ky', 'f', plot_format)
        singlePlot2D(kygrid, fgrid, dens_fky, 'dens', title, filename, 'ky', 'f', plot_format)

