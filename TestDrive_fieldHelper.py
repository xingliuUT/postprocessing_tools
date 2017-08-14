#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
from fieldHelper import *
from parIOHelper import *
import sys

suffix = sys.argv[1]

if not suffix =='.dat':
   suffix = '_'+suffix

pars = init_read_parameters(suffix)
nz = pars['nz0']
nx = pars['nx0']
show_plots = True
phi, apar = global_eigenfunctions(pars,suffix,show_plots,setTime=-1)

