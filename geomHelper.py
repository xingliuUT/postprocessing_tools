from parIOHelper import init_read_parameters
from read_write_geometry import *
import matplotlib.pyplot as plt
import numpy as np

def init_read_geometry(suffix, pars):

    geom_type = pars['magn_geometry'][1:-1]
    geom_file = geom_type + suffix
    geom_pars, geom_coeff = read_geometry_local(geom_file)

    return geom_type, geom_pars, geom_coeff

def local_grid_points(geom_coeff, show_plots = False):

    Z = geom_coeff['gl_z']
    R = geom_coeff['gl_R']

    if show_plots:
        plt.scatter(R,Z)
        plt.xlabel('R')
        plt.ylabel('Z')
        plt.axis('equal')
        plt.title('simulation grid points')
        plt.show()

    return R, Z

def ky(pars, geom_coeff, ktheta_cm = 1., show_plots = False):

    ggxx = geom_coeff['ggxx']
    ggxy = geom_coeff['ggxy']
    ggyy = geom_coeff['ggyy']

    gamma1 = ggxx * ggyy - ggxy ** 2

    kymin = pars['kymin']
    ky = np.sqrt(gamma1/ggxx)*kymin
    ky /= ky[pars['nz0']/2]
    ky *= ktheta_cm * 100.


    R, Z = local_grid_points(geom_coeff, show_plots)

    if show_plots:
        plt.plot(ky,label='ky (1/m)')
        plt.legend()
        plt.show()

    return ky


