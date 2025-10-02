import numpy as np
import h5py as h5py

def reset_snapshot_temperatures(filename=None):

    gamma_eos = 5./3.
    T_gas = 1.0e4       # K
    T_halo = 1.0e6      # K
    boltzmann_cgs = 1.380658e-16 # erg K^{-1}
    m_H = 1.6733e-24    # g
    mu = 0.62111795     # dimensionless mean molecular weight
    mean_molecular_weight = mu * m_H  # fully ionized solar metallicity gas
    c_v = (boltzmann_cgs / mean_molecular_weight) * (1.0 / (gamma_eos - 1.0))

    # conversion factor from specific energy in code units to specific energy in cgs
    unittime_cgs = 3.08568e16    # s (1 Gyr in s)
    unitlength_cgs = 3.085678e21 # cm (1 kpc in cm)
    unitenergypermass_cgs = unitlength_cgs**(2) * unittime_cgs**(-2) # erg g^{-1}

    # gas temperature if R < 20 kpc and |z| < 3 kpc
    u_diskgas = c_v * T_gas / unitenergypermass_cgs    # (specific internal energy)
    # gas temperature (otherwise)
    u_halogas = c_v * T_halo / unitenergypermass_cgs   # (specific internal energy)

    #%---- System of units used in GIZMO parameterfile
    #UnitLength_in_cm            3.085678e21    % 1.0 kpc
    #UnitMass_in_g               1.989e43  	    % 1.0e10 solar masses
    #UnitVelocity_in_cm_per_s    1.0e5   	    % 1 km/sec
    #UnitMagneticField_in_gauss  1.0   	        % 1 gauss

    # open HDF5 file for writing
    file = h5py.File(filename,'r+') 
    gas = file['PartType0/Coordinates']

    # Gas particles (particle type 0)
    x_g=gas[:,0]
    y_g=gas[:,1]
    z_g=gas[:,2]

    # set the initial internal energy per unit mass
    R_g = np.sqrt(x_g**2 + y_g**2)
    is_disk_gas = np.logical_and( R_g <= 20.0, np.abs(z_g) <= 3.0 )
    u_g = np.zeros(len(gas))
    u_g[is_disk_gas] = u_diskgas
    u_g[~is_disk_gas] = u_halogas

    # save the new specific internal energies
    file['PartType0/InternalEnergy'][:] = u_g

    # close the HDF5 file, which saves these outputs
    file.close()


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("filename")
    args = parser.parse_args()
    
    reset_snapshot_temperatures(args.filename)
