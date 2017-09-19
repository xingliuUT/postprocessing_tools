import matplotlib.pyplot as plt
import numpy as np
from geomHelper import *
from parIOHelper import *
import sys

suffix = sys.argv[1]

if not suffix =='.dat':
   suffix = '_'+suffix

pars = init_read_parameters(suffix)

geom_type, geom_pars, geom_coeff = init_read_geometry(suffix, pars)

zGrid(geom_coeff, pars)
