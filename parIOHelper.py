from ParIO import *

def init_read_parameters_file(suffix):

    par = Parameters()
    par.Read_Pars('parameters'+suffix)
    pars = par.pardict

    return pars
