import argparse
from matplotlib.cm import colors
import numpy as np
from gizmo_analysis import *

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='plot a GIZMO output')
    parser.add_argument('input_filename', help='hdf5 snapshot filename')
    parser.add_argument('output_filename', help='output plot file')
    args = parser.parse_args()

    snap = args.input_filename
    output_file = args.output_filename

    pdata = load_hydro_data(snap)
    temp = compute_temperature(pdata)
    mesh = compute_mesh(pdata)

    import matplotlib.pyplot as plt

    ## make phase plot! (temperature vs n_H)
    dens = mesh.Density() * unitdensity_per_H

    lognHmin = 0.0
    lognHmax = 5.0
    logTmin = 0.0
    logTmax = 5.0
    
    h = plt.hist2d(np.log10(dens), np.log10(temp),
                    bins=100, norm=colors.LogNorm(), density=True,
                    range = [[lognHmin,lognHmax], [logTmin,logTmax]])

    plt.colorbar(label=r'probability density by mass')
    plt.xlabel(r'$\log_{10}$ density (H cm$^{-3}$)')
    plt.ylabel(r'$\log_{10}$ temperature (K)')
    plt.savefig(output_file, dpi=fig_dpi)
    plt.close()