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
field = fieldfile('field'+suffix,pars)
#field_step_time(field)
momen = momfile('mom_e'+suffix,pars)
#momen_step_time(momen)

nz = pars['nz0']
nx = pars['nx0']
dz = 2.0/nz
zgrid = np.arange(nz)/float(nz-1)*(2.0-dz)-1.0
if 'lx_a' in pars:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx_a']+pars['x0']-pars['lx_a']/2.0
else:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx'] - pars['lx']/2.0

show_plots = True
plot_format = 'display'

kygrid = range(0, pars['nky0'])
if 1 == 1:
    time, phi, apar = reflectometer(field, nz / 2, kygrid, -1)
    if show_plots:
        title = 'time =' + str(np.round(time,1)) + ' all ky\'s'
        filename = 'reflectometer_phi_Apar.ps'
        doublePlot1D(xgrid, phi, apar, 'phi', 'apar', title, filename, plot_format)
#ky = 5

#for ky in kygrid:
if 1 == 0:
    time, phi, apar = global_eigenfunctions(field, nz / 2, ky, -1, setTime=-1)
    if show_plots:
        title = 'time =' + str(np.round(time,1)) + ', ky index =' + str(ky)
        filename = 'phi_Apar_ky='+str(ky)+'.ps'
        doublePlot1D(xgrid, phi, apar, 'phi', 'apar', title, filename, plot_format)
if 1 == 0:
    tStart = 350.
    tEnd = 390.
    phi, apar = t_avg_global_eigenfns(field,ky,tStart,tEnd)
    if show_plots:
        title = 'time=' + str(tStart) + '~' + str(tEnd) + ', ky index =' + str(ky)
        filename = 't_avg_phi_Apar_ky='+str(ky)+'.ps'
        doublePlot2D(xgrid, zgrid, phi, apar, 'phi', 'apar', title, filename, plot_format)
if 1 == 0:
    tStart = 350.
    tEnd = 390.
    dens, tperp = t_avg_global_moms(momen,ky,tStart,tEnd)
    if show_plots:
        title = 'time=' + str(tStart) + '~' + str(tEnd) + ', ky index =' + str(ky)
        filename = 't_avg_dens_Tperp_ky='+str(ky)+'.ps'
        doublePlot2D(xgrid, zgrid, dens, tperp, 'dens', 'tperp', title, filename, plot_format)

if 1 == 0:
#kygrid = range(0, pars['nky0'])
#for ky in kygrid:
    time, phi, apar = global_eigenfunctions(field,ky,setTime=-1)
    if show_plots:
        title = 'time =' + str(np.round(time,1)) + ', ky index =' + str(ky)
        filename = 'phi_Apar_ky='+str(ky)+'.ps'
        doublePlot2D(xgrid, zgrid, phi, apar, 'phi', 'apar', title, filename, 'ps')
    time, dens, tperp = global_moments(momen,ky,setTime=-1)
    if show_plots:
        title = 'time =' + str(np.round(time,1)) + ', ky index =' + str(ky)
        filename = 'dens_Tperp_ky='+str(ky)+'.ps'
        doublePlot2D(xgrid, zgrid, dens, tperp, 'n', 'Tperp', title, filename, 'ps')

