import numpy as np
import argparse
from pathlib import Path
from gizmo_analysis import *


def if_not_exists(filename, func):
    if not Path(filename).exists():
        func()

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='plot a GIZMO output')
    parser.add_argument('filename', help='hdf5 snapshot filename')
    args = parser.parse_args()
    
    snaps = [Path(args.filename)]    
    for snap in snaps:
        
        pdata = load_hydro_data(snap)
        temp = compute_temperature(pdata)
        mesh = compute_mesh(pdata)
        stars = load_stars(snap)
        
        phaseplot_file = snap.parent.with_name(snap.name).with_suffix('.phaseplot.png')
        if_not_exists(phaseplot_file,
                      lambda: save_phase_plot(mesh.Density(), temp, phaseplot_file))

        rmax = 15.0 # kpc

        sliceplot_file = snap.parent.with_name(snap.name).with_suffix('.slice_temperature+stars.png')
        if_not_exists(sliceplot_file,
                      lambda: save_slice_plot(mesh, temp, sliceplot_file,
                            colorbar_label=r'Temperature (K)', rmax=rmax, star_coords=stars))

        projplot_file = snap.parent.with_name(snap.name).with_suffix('.projection_density+stars.png')
        if_not_exists(projplot_file,
                      lambda: save_density_projection_plot(mesh, projplot_file,
                                                           rmax=rmax,
                                                           star_coords=stars,
                                                           bfield=pdata['MagneticField']))

        sliceplot_file = snap.parent.with_name(snap.name).with_suffix('.slice_temperature.png')
        if_not_exists(sliceplot_file,
                      lambda: save_slice_plot(mesh, temp, sliceplot_file,
                            colorbar_label=r'Temperature (K)', rmax=rmax))

        projplot_file = snap.parent.with_name(snap.name).with_suffix('.projection_density.png')
        if_not_exists(projplot_file,
                      lambda: save_density_projection_plot(mesh, projplot_file,
                                                           rmax=rmax,
                                                           bfield=pdata['MagneticField']))
                      
