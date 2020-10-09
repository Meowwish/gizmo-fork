import numpy as np
import argparse
from pathlib import Path
from gizmo_analysis import *

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='plot a GIZMO output')
    parser.add_argument('filename', help='hdf5 snapshot filename')
    args = parser.parse_args()
    
    snaps = [Path(args.filename)]    
    for snap in snaps:

        phaseplot_file = snap.parent.with_name(snap.name).with_suffix('.phaseplot.png')
        sliceplot_file = snap.parent.with_name(snap.name).with_suffix('.slice_temperature.png')
        zsliceplot_file = snap.parent.with_name(snap.name).with_suffix('.zslice_density.png')
        ztempplot_file = snap.parent.with_name(snap.name).with_suffix('.zslice_temperature.png')
        projplot_file = snap.parent.with_name(snap.name).with_suffix('.projection_density.png')
        magpressure_sliceplot = snap.parent.with_name(snap.name).with_suffix('.magnetic_energy.png')

        files = [phaseplot_file, sliceplot_file, zsliceplot_file, projplot_file, ztempplot_file,
                 magpressure_sliceplot]
        files_exist = [f.exists() for f in files]

        if all(files_exist):
            continue # skip this snapshot

        ## plotting parameters
        rmax = 15.0 # kpc
        
        rawdata = load_hydro_data(snap)
        rawtemp = compute_temperature(rawdata)

        pdata, temp = apply_radius_cut(rawdata, rawtemp, rmax=np.sqrt(2.)*rmax)
        mesh = compute_mesh(pdata)
        #stars = load_stars(snap)

        if not phaseplot_file.exists() and rawtemp is not None:
            save_phase_plot(mesh.Density(), temp, phaseplot_file)

        if not sliceplot_file.exists() and rawtemp is not None:
            save_slice_plot(mesh, temp, sliceplot_file,
                            colorbar_label=r'Temperature (K)', rmax=rmax)


        if not magpressure_sliceplot.exists():
            bfield = pdata['MagneticField'] * unitbfield_cgs # gauss
            mag_energy_density = np.sqrt(np.einsum('ij,ij->i', bfield, bfield)) / (8.0*np.pi)
            print(f"mean magnetic energy density = {np.mean(mag_energy_density))}")
            save_slice_plot(mesh, mag_energy_density,
                            magpressure_sliceplot,
                            colorbar_label=r'magnetic energy density (erg cm$^{-3}$)',
                            rmax=rmax,
                            bfield=pdata['MagneticField'],
                            vmin=1.0e-9,
                            vmax=1.0e-15)

        if not zsliceplot_file.exists():
            save_slice_plot(mesh, mesh.Density()*unitdensity_per_H, zsliceplot_file,
                            colorbar_label=r'Density (g cm$^-3$)', rmax=rmax, plane='y',
                            vmin=1e-5, vmax=1e5)

        if not ztempplot_file.exists() and rawtemp is not None:
            save_slice_plot(mesh, temp, ztempplot_file,
                            colorbar_label=r'Temperature (K)', rmax=rmax, plane='y',
                            vmin=10., vmax=1.0e7)

        if not projplot_file.exists():
            save_density_projection_plot(mesh, projplot_file,
                                         rmax=rmax,
                                         bfield=pdata['MagneticField'])
        
