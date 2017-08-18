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
#plot_format = 'display'
plot_format = 'ps'
nf = 200
lf = 10.

tStart = 350.
tEnd = 390.
kygrid = range(pars['nky0'])
#kygrid = [0]
zInd = nz/2
kyInd = -1
xInd = nx*5/8
ky = 0
if 1 == 0:
    time, phi, apar = global_eigenfunctions(field, -1, ky, -1)
    if show_plots:
        title = 'time =' + str(np.round(time,1)) + ', ky index =' + str(ky)
        filename = 'phi_Apar_ky='+str(ky)+'.ps'
        doublePlot2D(xgrid, zgrid, phi, apar, 'phi', 'apar', title, filename, 'x', 'z',plot_format)
if 1 == 1:
    #time, dens, tperp = global_moments(momen,zInd,kyInd,xInd)
    #if show_plots:
    if 1 == 0:
        title = 'time =' + str(np.round(time,1))
        filename = 'dens_Tperp_.ps'
        doublePlot1D(kygrid, dens, tperp, 'n', 'Tperp', title, filename, plot_format)

    tgrid, dens_tky, tperp_tky = momen_tky(momen,zInd,kyInd,xInd,tStart,tEnd)
    if 1 == 0:
        title = 'xgrid: time, ygrid: ky, x ='+str(xgrid[xInd])
        filename = 'tky_dens_tperp.ps'
        doublePlot2D(kygrid, tgrid, dens_tky, tperp_tky, 'dens', 'tperp', title, filename, 'ky', 't', plot_format)

    fgrid, dens_fky = momen_fx(dens_tky, tgrid, nf, lf)
    fgrid, tperp_fky = momen_fx(tperp_tky, tgrid, nf, lf)
    if 1 == 1:
        title = 'xgrid: frequency, ygrid: x, all ky\'s added'
        filename = 'fky_dens_tperp02.ps'
        doublePlot2D(kygrid, fgrid, dens_fky, tperp_fky, 'dens', 'tperp', title, filename, 'ky', 'f', plot_format)

if 1 == 0:
    tgrid, dens, tperp = momen_tx(momen, nz / 2, kygrid, -1, tStart, tEnd)
#    for i in range(nx/2, nx, 16):
    if 1 == 0:
        f_grid, dens_f = windowFFT(np.array(tgrid), dens[:,i], nf, lf, 'dens, x ='+str(xgrid[i]))
        f_grid, tperp_f = windowFFT(np.array(tgrid), tperp[:,i], nf, lf, 'Tperp, x ='+str(xgrid[i]))
    if 1 == 1:
        title = 'xgrid: time, ygrid: x, all ky\'s added'
        filename = 'tx_dens_tperp.ps'
        doublePlot2D(xgrid, tgrid, dens, tperp, 'dens', 'tperp', title, filename, 'rhot', 't', plot_format)
        fgrid, dens_fx = momen_fx(dens, tgrid, nf, lf)
        fgrid, tperp_fx = momen_fx(tperp, tgrid, nf, lf)
        title = 'xgrid: frequency, ygrid: x, all ky\'s added'
        filename = 'fx_dens_tperp.ps'
        doublePlot2D(xgrid, fgrid, dens_fx, tperp_fx, 'dens', 'tperp', title, filename, 'rhot', 't', plot_format)

if 1 == 0:
    tgrid, phi, apar = field_tx(field, nz / 2, kygrid, -1, tStart, tEnd)
    for i in range(nx):
#    if 1 == 0:
        windowFFT(np.array(tgrid), phi[:,i], 'phi, x ='+str(xgrid[i]))
        windowFFT(np.array(tgrid), apar[:,i], 'apar, x ='+str(xgrid[i]))
    if 1 == 1:
        title = 'xgrid: time, ygrid: x, all ky\'s added'
        filename = 'tx_phi_Apar.ps'
        #for i in range(len(tgrid)):
        #    title = 'all ky\'s added, time =' + str(tgrid[i])
        #    doublePlot1D(xgrid, phi[i,:], apar[i,:], 'phi', 'apar', title, filename, plot_format)
        doublePlot2D(xgrid, tgrid, phi, apar, 'phi', 'apar', title, filename, 'rhot', 't', plot_format)
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

