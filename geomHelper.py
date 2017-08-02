from parIOHelper import init_read_parameters
from read_write_geometry import *

def init_read_geometry(suffix, pars):

    geom_type = pars['magn_geometry'][1:-1]
    geom_file = geom_type + suffix
    geom_pars, geom_coeff = read_geometry_local(geom_file)

    return geom_type, geom_pars, geom_coeff

