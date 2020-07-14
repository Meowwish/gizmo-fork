import numpy as np
import yt
from pathlib import Path

if __name__ == '__main__':

    output_dir = Path("./output")
    snaps = list(output_dir.glob("snapshot_*.hdf5"))
    bbox = 1000.
    bounds = [[-bbox,bbox],[-bbox,bbox],[-bbox,bbox]]

    for snap in snaps:
        ds = yt.load(str(snap), bounding_box=bounds)
        yt.ProjectionPlot(ds, 'z', 'density').set_width((20,'kpc')).save()
    