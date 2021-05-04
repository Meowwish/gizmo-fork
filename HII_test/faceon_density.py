import argparse
import matplotlib
import scipy.interpolate
import scipy.optimize
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

    # compute estimate of r_HII
    Lbox = 20.
    lbox_cgs = Lbox * unitlength_cgs
    massbox_cgs = np.sum(pdata['Masses']) * unitmass_cgs
    volbox_cgs = lbox_cgs**3
    mean_density_cgs = massbox_cgs / volbox_cgs
    n_H = hydrogen_massfrac * mean_density_cgs / m_H
    N_photons = 1.0e49 # s^{-1}
    beta = 3.0e-13
    r1_approx = (3.0 * N_photons / (4.0 * np.pi * n_H * n_H * beta))**(1./3.)
    print(f"r1_approx: {(r1_approx / unitlength_cgs):.4f} pc")


    Tfin = 1.0e4 # Kelvins
    molw_i = 4.0 / (8 - 5 * (1 - hydrogen_massfrac))
    e_ion = boltzmann_cgs * Tfin / ((gamma - 1) * molw_i * m_H) / unitvelocity_cgs**2
    cs = np.sqrt(boltzmann_cgs * Tfin / (molw_i*m_H)) / 1.0e5
    print(f"sound speed: {cs} km/s")

    eps = 5.0e-2
    HII_region_mask = np.abs((pdata['InternalEnergy'] - e_ion)/e_ion) < eps
    HII_region_field = np.zeros_like(temp)
    HII_region_field[HII_region_mask] = 1.0     # select HII regions

    print(f"e_ion = {e_ion:.3e}")
    print(f"HII particles: {np.count_nonzero(HII_region_mask)}")

    x,y,field_slice = save_slice_plot(mesh, HII_region_field, args.output_filename,
                    colorbar_label=r"HII region",
                    bfield=None, star_coords=None, xymin=0., xymax=20., plane='z',
                    norm=matplotlib.colors.Normalize(vmin=0.0, vmax=1.0))

    #import matplotlib.pyplot as plt
    #ax = plt.gca()
    #circle = matplotlib.patches.Circle((Lbox/2.,Lbox/2.), radius=r1_approx,
    #            facecolor='none', edgecolor='white', linewidth=3, zorder=100)
    #ax.add_patch(circle)

    # RectBivariateSpline
    interpolant = scipy.interpolate.RectBivariateSpline(x,y,field_slice)
    midpoint = 10.

    f = lambda y: interpolant(midpoint, y) - 0.5
    r_HII = scipy.optimize.bisect(f, midpoint, 20.) - midpoint
    print(f"radius of HII region: {r_HII:.4f} pc")