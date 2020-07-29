import numpy as np
import yt
from pathlib import Path

if __name__ == '__main__':

    output_dir = Path("./output")
    snaps = list(output_dir.glob("snapshot_*.hdf5"))
    bbox = 1000.
    bounds = [[-bbox,bbox],[-bbox,bbox],[-bbox,bbox]]

    for snap in snaps:
        #ds = yt.load(str(snap), bounding_box=bounds, kernel_name='cubic', over_refine_factor=3, n_ref=16)
        ds = yt.load(str(snap), bounding_box=bounds, kernel_name='cubic', over_refine_factor=2)
        #ds = yt.load(str(snap), bounding_box=bounds, kernel_name='cubic')

        plt = yt.ProjectionPlot(ds, 'z', 'density')
        plt.set_width((40,'kpc'))
        if ('PartType4', 'StellarFormationTime') in ds.field_list:
            plt.annotate_particles(1.0, ptype='PartType4', p_size=5.0)
            plt.annotate_text((-17.5,17.5), 'Projected gas density + stars overplotted', coord_system='plot')

        plt.save()

        plt = yt.SlicePlot(ds, 'z', 'temperature')
        plt.set_width((40,'kpc'))
        plt.save()

        my_sphere = ds.sphere("c", (50, "kpc"))
        plot = yt.PhasePlot(my_sphere, "density", "temperature", ["cell_mass"],
                            weight_field=None)
        plot.set_unit('density', 'Msun/pc**3')
        plot.set_unit('cell_mass', 'Msun')
        #plot.set_ylim(1e-7, 1e7)
        plot.save()
        
        #yt.SlicePlot(ds, 'x', 'temperature').set_width((40,'kpc')).save()
    
