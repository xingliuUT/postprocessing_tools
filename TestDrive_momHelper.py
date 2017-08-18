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

pars = init_read_parameters(suffix)
field = fieldfile('field'+suffix,pars)
#field_step_time(field)
momen = momfile('mom_e'+suffix,pars)
#momen_step_time(momen)
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
plot_format = 'display'
#plot_format = 'ps'
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

geom_type, geom_pars, geom_coeff = init_global_geometry(suffix, pars)
q, Cy = q_Cy(geom_coeff)
qCy = np.array(q * Cy)
#print np.shape(q), np.shape(Cy), np.shape(qCy)
#print np.shape(zgrid)
ymatrix = np.outer(zgrid*np.pi, qCy)
#print np.shape(ymatrix)
#filename = '0'
#doublePlot1D(xgrid, q, Cy, 'q', 'Cy', 'geometry coeff', filename, plot_format)
#if 1 == 1:
dens_xz = np.zeros((nz, nx),dtype='complex128')
dens_naive = np.zeros((nz, nx),dtype='complex128')
for ky in kygrid:
    time, this_dens, this_tperp = global_moments(momen, -1, ky, -1)
    dens_xz += np.multiply(this_dens, np.exp(zi * ky * ymatrix))
    dens_naive += this_dens
#    print np.shape(this_dens)
#    if show_plots:
if 1 == 1:
    title = 'time =' + str(np.round(time,1))
    filename = 'dens01.ps'
    doublePlot2D(xgrid, zgrid, dens_xz, dens_naive, 'dens_xz', 'dens_naive', title, filename, 'x', 'z',plot_format)
