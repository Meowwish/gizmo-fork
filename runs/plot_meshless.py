import numpy as np
from pathlib import Path
from gizmo_analysis import *

def if_not_exists(filename, func):
    if not Path(filename).exists():
        func()

if __name__ == '__main__':

    output_dir = Path("./output")
    snaps = list(output_dir.glob("snapshot_102.hdf5"))

    for snap in snaps:
        print(f"plotting {snap}...")
        pdata = load_hydro_data(snap)
        temp = compute_temperature(pdata)
        mesh = compute_mesh(pdata)
        stars = load_stars(snap)
        if stars is not None:
            print(f"\tNumber of stars = {len(stars)}")
        
        phaseplot_file = snap.parent.with_name(snap.name).with_suffix('.phaseplot.png')
        if_not_exists(phaseplot_file,
                      lambda: save_phase_plot(mesh.Density(), temp, phaseplot_file))

        rmax = 20.0 # kpc

        sliceplot_file = snap.parent.with_name(snap.name).with_suffix('.slice_temperature+stars.png')
        if_not_exists(sliceplot_file,
                      lambda: save_slice_plot(mesh, temp, sliceplot_file,
                            colorbar_label=r'Temperature (K)', rmax=rmax, star_coords=stars))

        projplot_file = snap.parent.with_name(snap.name).with_suffix('.projection_density+stars.png')
        if_not_exists(projplot_file,
                      lambda: save_density_projection_plot(mesh, projplot_file,
                                                           rmax=rmax, star_coords=stars))

        sliceplot_file = snap.parent.with_name(snap.name).with_suffix('.slice_temperature.png')
        if_not_exists(sliceplot_file,
                      lambda: save_slice_plot(mesh, temp, sliceplot_file,
                            colorbar_label=r'Temperature (K)', rmax=rmax))

        projplot_file = snap.parent.with_name(snap.name).with_suffix('.projection_density.png')
        if_not_exists(projplot_file,
                      lambda: save_density_projection_plot(mesh, projplot_file,
                                                           rmax=rmax))
                      
