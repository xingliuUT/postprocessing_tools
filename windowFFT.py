import numpy as np
import matplotlib.pyplot as plt

zi = complex(0, 1)

def windowFFT(tgrid, field, nf, lf, tt, show_plots = True):
    nf = int(nf)
    lf = float(lf)
    tgrid_n = (tgrid - tgrid[0]) / (tgrid[-1] - tgrid[0])
    window = np.cos(np.pi * tgrid_n - np.pi / 2.)
    if 1 == 0:
        plt.plot(tgrid, abs(field)/max(abs(field)))
        plt.plot(tgrid, window)
        plt.show()

    fgrid = np.linspace(-lf,lf,nf,endpoint = False)
    field_f = np.empty(0, dtype = 'complex128')
    for f in fgrid:
        this_field_f = 0.
        for i in range(len(field) - 1):
            this_field_f += 0.5*(field[i]*np.exp(zi*f*tgrid[i]) + \
                          field[i + 1]*np.exp(zi*f*tgrid[i + 1])) * \
                          (tgrid[i + 1] - tgrid[i])
        field_f = np.append(field_f,this_field_f)
    if show_plots:
        plt.plot(fgrid, np.real(field_f), '+-')
        plt.plot(fgrid, np.imag(field_f), '+-')
        plt.title(tt)
        plt.show()

