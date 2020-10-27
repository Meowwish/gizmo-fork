import numpy as np
import argparse
from pathlib import Path
from gizmo_analysis import *
import yt
import yt.visualization.api
import yt.units

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='plot a GIZMO output')
    parser.add_argument('filename', help='hdf5 snapshot filename')
    args = parser.parse_args()
    
    snaps = [Path(args.filename)]    
    for snap in snaps:

        phaseplot_file = snap.with_suffix('.phaseplot.png')
        sliceplot_file = snap.with_suffix('.slice_temperature.png')
        zsliceplot_file = snap.with_suffix('.zslice_density.png')
        ztempplot_file = snap.with_suffix('.zslice_temperature.png')
        projplot_file = snap.with_suffix('.projection_density.png')
        zprojplot_file = snap.with_suffix('.zprojection_density.png')
        magpressure_sliceplot = snap.with_suffix('.magnetic_energy.png')

        files = [projplot_file, zprojplot_file]
        files_exist = [f.exists() for f in files]

        if all(files_exist):
            continue # skip this snapshot

        ## plotting parameters
        rmax = 15.0 # kpc
        plot_zslice = True
        
        rawdata = load_hydro_data(snap)
        rawtemp = compute_temperature(rawdata)
        pdata, temp = apply_radius_cut(rawdata, rawtemp, rmax=np.sqrt(2.)*rmax)
        mesh = compute_mesh(pdata)

        #if not phaseplot_file.exists() and rawtemp is not None:
        #    save_phase_plot(mesh.Density(), temp, phaseplot_file)

        if not projplot_file.exists():
            save_density_projection_plot(mesh, projplot_file,
                                         rmax=rmax,
                                         bfield=None,
                                         projection_length=1.0)
        
        if not zprojplot_file.exists():
            save_zdensity_projection_plot(mesh, zprojplot_file,
                                         rmax=rmax,
                                         bfield=None,
                                         plane='y',
                                         projection_length=0.5,
                                         vmin=1e-4, vmax=1e2)

        # load data with yt
        # bbox_lim = 1e5 #kpc
        # bbox = [[-bbox_lim,bbox_lim],
        #         [-bbox_lim,bbox_lim],
        #         [-bbox_lim,bbox_lim]]
        # ds = yt.load(str(snap), bounding_box=bbox)

        # ds.add_field('Bx', function=lambda field, data: data[('PartType0','MagneticField')][:,0]) 
        # ds.add_field('By', function=lambda field, data: data[('PartType0','MagneticField')][:,1]) 
        # ds.add_field('Bz', function=lambda field, data: data[('PartType0','MagneticField')][:,2]) 

        # num_streamlines = 10
        # scale = 10.0*yt.units.kpc
        # start_dx = np.random.random((num_streamlines, 3))*scale - 0.5*scale
        # start = ds.domain_center + start_dx

        # streamlines = yt.visualization.api.Streamlines(ds, start, 'Bx', 'By', 'Bz',
        #                                                length=1.0*yt.units.kpc)
        # print("integrating streamlines...")
        # streamlines.integrate_through_volume()
        # lines = streamlines.streamlines
        # print("done.")

        ## the first-created slice plot is very computationally expensive (why??)

        # if not magpressure_sliceplot.exists():
        #     bfield = pdata['MagneticField'] * unitbfield_cgs # gauss
        #     mag_energy_density = np.sqrt(np.einsum('ij,ij->i', bfield, bfield)) / (8.0*np.pi)
        #     save_slice_plot(mesh, mag_energy_density,
        #                     magpressure_sliceplot,
        #                     colorbar_label=r'magnetic energy density (erg cm$^{-3}$)',
        #                     rmax=rmax,
        #                     bfield=pdata['MagneticField'],
        #                     vmin=1.0e-9,
        #                     vmax=1.0e-15)

        #if not sliceplot_file.exists() and rawtemp is not None:
        #    save_slice_plot(mesh, temp, sliceplot_file,
        #                    colorbar_label=r'Temperature (K)', rmax=rmax)
        
        # if not zsliceplot_file.exists() and plot_zslice:
        #     save_slice_plot(mesh, mesh.Density()*unitdensity_per_H, zsliceplot_file,
        #                     colorbar_label=r'Density (g cm$^-3$)', rmax=rmax, plane='y',
        #                     vmin=1e-5, vmax=1e5, bfield=pdata['MagneticField'])

        # if not ztempplot_file.exists() and rawtemp is not None and plot_zslice:
        #     save_slice_plot(mesh, temp, ztempplot_file,
        #                     colorbar_label=r'Temperature (K)', rmax=rmax, plane='y',
        #                     vmin=10., vmax=1.0e7, bfield=pdata['MagneticField'])

