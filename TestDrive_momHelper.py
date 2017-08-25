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
geom_type, geom_pars, geom_coeff = init_global_geometry(suffix, pars)

zi = complex(0, 1)
nz = pars['nz0']
nx = pars['nx0']
dz = 2.0/nz
zgrid = np.arange(nz)/float(nz-1)*(2.0-dz)-1.0
if 'lx_a' in pars:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx_a']+pars['x0']-pars['lx_a']/2.0
else:
    xgrid = np.arange(nx)/float(nx-1)*pars['lx'] - pars['lx']/2.0

show_plots = True
show_xz = False
plot_format = 'display'
#plot_format = 'ps'
nf = 200
lf = 10.

kygrid = range(pars['nky0'])
zInd = nz/2
kyInd = -1

if 1 == 1:
    tgrid, dens_tx, tperp_tx = momen_tx(momen, \
                  geom_coeff, \
                  zgrid, \
                  kygrid, \
                  xgrid, \
                  zInd, \
                  tStart, \
                  tEnd, \
                  show_xz, \
                  plot_format)
#    np.savetxt('dens_tx.txt', dens_tx.view(float))
#    np.savetxt('tperp_tx.txt', tperp_tx.view(float))
#    np.savetxt('tgrid.txt', tgrid)
    if show_plots:
#        plot_format = 'ps'
        title = ' '
        filename = 'tx_dens_tperp.ps'
        doublePlot2D(xgrid, tgrid, dens_tx, tperp_tx, 'dens_tx', 'tperp_tx', title, filename, 'x', 't',plot_format)
if 1 == 1:
    fgrid, dens_fx = momen_fx(dens_tx, tgrid, nf, lf)
    fgrid, tperp_fx = momen_fx(tperp_tx, tgrid, nf, lf)
#    title = 'tStart='+str(tStart)+', tEnd='+str(tEnd)
    np.savetxt('dens_fx.txt', dens_fx.view(float))
    np.savetxt('tperp_fx.txt', tperp_fx.view(float))
    np.savetxt('fgrid.txt', fgrid)
    title = ' '
    filename = 'phi_fx01.ps'
#    singlePlot2D(xgrid, fgrid, dens_fx, 'phi_fx', title, filename, 'x', 'f',plot_format)
    doublePlot2D(xgrid, fgrid, dens_fx, tperp_fx, 'phi_fx', 'apar_fx', title, filename, 'x', 'f',plot_format)

    for i in range(nx / 2, nx * 3 / 4, 8):
        dens_f = dens_fx[:,i]

        plt.plot(fgrid, abs(dens_f))
        plt.ylabel('abs phi')
        plt.xlabel('f')
        plt.title(str(i))
        plt.show()

