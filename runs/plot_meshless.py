import numpy as np
from pathlib import Path
from gizmo_analysis import *

if __name__ == '__main__':

    output_dir = Path("./output")
    snaps = list(output_dir.glob("snapshot_*.hdf5"))

    for snap in snaps:
        pdata = load_hydro_data(snap)
        temp = compute_temperature(pdata)
        mesh = compute_mesh(pdata)
        
        save_phase_plot(mesh.Density(), temp,
                        snap.parent.with_name(snap.name).with_suffix('.phaseplot.pdf'))
        save_slice_plot(mesh, temp,
                        snap.parent.with_name(snap.name).with_suffix('.slice_temperature.pdf'),
                        colorbar_label=r'Temperature (K)')
        save_density_projection_plot(mesh,
                        snap.parent.with_name(snap.name).with_suffix('.projection_density.pdf'))