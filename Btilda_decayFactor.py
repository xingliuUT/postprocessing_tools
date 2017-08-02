import sys
from parIOHelper import *
from geomHelper import *

suffix = sys.argv[1]
if not suffix == '.dat':
    suffix = '_' + suffix

efit_file_name = sys.argv[2]
ktheta_cm = 0.65
calcRZ = False # calculate Mirnov coil pos. based on exp.

pars = init_read_parameters(suffix)

geom_type, geom_pars, geom_coeff = init_read_geometry(suffix, pars)

R, Z = local_grid_points(geom_coeff, True)

ky_fluxsurface = ky(pars, geom_coeff, ktheta_cm, False)

if calcRZ:    # if this is the local linear run at rhot = 0.999...
    topInd = np.argmax(Z)
    botInd = np.argmin(Z)
    Z1 = Z[botInd : topInd + 1]
    mirnovZ = 0.1
    zInd = np.argmin(abs(Z1 - mirnovZ)) + botInd
    mirnovZ = Z[zInd]
    mirnovR = R[zInd]
else:    # if the probe's position at LCFS is found
    mirnovR = 0.86669 + 0.02
    mirnovZ = 0.1015
print('Mirnov probe is at ({:.4f}, {:.4f})'.format(mirnovR, mirnovZ))

print suffix, efit_file_name
