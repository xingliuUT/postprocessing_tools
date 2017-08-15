#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from fieldHelper import *
from momHelper import *
from parIOHelper import *
from plotHelper import *
import sys

suffix = sys.argv[1]

if not suffix =='.dat':
   suffix = '_'+suffix

pars = init_read_parameters(suffix)

nz = pars['nz0']
nx = pars['nx0']
dz = 2.0/nz
zgrid = np.arange(nz)/float(nz-1)*(2.0-dz)-1.0
if 'lx_a' in pars:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx_a']+pars['x0']-pars['lx_a']/2.0
else:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx'] - pars['lx']/2.0

show_plots = True
#kygrid = range(0, pars['nky0'])
kygrid = [5, 6]
for ky in kygrid:
    time, phi, apar = global_eigenfunctions(pars,suffix,ky,show_plots,setTime=-1)
    if show_plots:
        title = 'time =' + str(np.round(time,1)) + ', ky index =' + str(ky)
        filename = 'phi_Apar_ky='+str(ky)+'.ps'
        doublePlot2D(xgrid, zgrid, phi, apar, 'phi', 'apar', title, filename, 'ps')
    time, dens, tperp = global_moments(pars,suffix,ky,show_plots,setTime=-1)
    if show_plots:
        title = 'time =' + str(np.round(time,1)) + ', ky index =' + str(ky)
        filename = 'dens_Tperp_ky='+str(ky)+'.ps'
        doublePlot2D(xgrid, zgrid, dens, tperp, 'n', 'Tperp', title, filename, 'ps')

