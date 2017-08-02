import sys
from parIOHelper import *
from geomHelper import *

suffix = sys.argv[1]
if not suffix == '.dat':
    suffix = '_' + suffix

efit_file_name = sys.argv[2]

pars = init_read_parameters(suffix)

geom_type, geom_pars, geom_coeff = init_read_geometry(suffix, pars)

ky_fluxsurface = ky(pars, geom_coeff, True)

print suffix, efit_file_name
