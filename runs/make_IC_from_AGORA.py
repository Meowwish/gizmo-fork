import numpy as np
import h5py as h5py
import pandas as pd # to read csv files
from pathlib import Path

def make_agora_IC(filename, magnetic_fields=True):

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


    # UNITS in csv files
    # ---------------
    # Velocity: km/s
    # Mass: 10^9 Msun
    # Length: kpc

    agora_dir = Path("../../agora_highres")
    
    # dark matter
    dm = pd.read_csv(agora_dir / "halo.dat", delimiter=' ', skipinitialspace=True,
                        names=['x','y','z','vx','vy','vz','mass'])
    # bulge stars
    bulge = pd.read_csv(agora_dir / "bulge.dat", delimiter=' ', skipinitialspace=True,
                        names=['x','y','z','vx','vy','vz','mass'])
    # disk stars
    disk = pd.read_csv(agora_dir / "disk.dat", delimiter=' ', skipinitialspace=True,
                        names=['x','y','z','vx','vy','vz','mass'])
    # gas particles
    gas = pd.read_csv(agora_dir / "gas.dat", delimiter=' ', skipinitialspace=True,
                        names=['x','y','z','vx','vy','vz','mass','u_gas'])

    #%---- System of units used in GIZMO parameterfile
    #UnitLength_in_cm            3.085678e21    % 1.0 kpc
    #UnitMass_in_g               1.989e43  	    % 1.0e10 solar masses
    #UnitVelocity_in_cm_per_s    1.0e5   	    % 1 km/sec
    #UnitMagneticField_in_gauss  1.0   	        % 1 gauss

    # WARNING: need to convert mass units from (10^9 Msun) -> (10^10 Msun)
    mass_fac = 1.0/10.0     # ( 10^9 / 10^10 )


    ## Gas particles (particle type 0)

    # position
    x_g=gas['x']
    y_g=gas['y']
    z_g=gas['z']
    R = np.sqrt( x_g**2 + y_g**2 ) # kpc
    z = z_g                     # kpc

    # velocity
    vx_g=gas['vx']
    vy_g=gas['vy']
    vz_g=gas['vz']

    # magnetic fields
    # (by default, identically zero)
    bx_g = 0.
    by_g = 0.
    bz_g = 0.

    if (magnetic_fields == True):
        print("Generating magnetic fields...")
        R_scale = 3.43218           # kpc
        z_scale = 0.343218          # kpc
        B0 = 10.e-6                 # gauss
        B_phi = B0 * np.exp(-R/R_scale) * np.exp(-abs(z)/z_scale)
        theta = np.arctan2(y_g, x_g)   
        # re-define magnetic fields
        bx_g = B_phi * (-1.0 * np.sin(theta))
        by_g = B_phi * (np.cos(theta))
        bz_g = 0.

    # masses
    m_g=gas['mass']*mass_fac
    print(f"Gas particle mass: {1e10*m_g[0]:.3g}")

    # specific internal gas energy [i.e., NOT including magnetic energy]
    R_g = np.sqrt(x_g**2 + y_g**2)
    is_disk_gas = np.logical_and( R_g <= 20.0, np.abs(z_g) <= 3.0 )
    u_g = np.zeros(len(gas))
    u_g[is_disk_gas] = u_diskgas
    u_g[~is_disk_gas] = u_halogas

    # particle ID
    Ngas = len(gas)
    print(f"Number of gas particles: {Ngas:.3g}")
    id_g=np.arange(1, Ngas+1)


    ## DM particles (particle type 1)

    # position
    x_dm = dm['x']
    y_dm = dm['y']
    z_dm = dm['z']

    # velocities
    vx_dm = dm['vx']
    vy_dm = dm['vy']
    vz_dm = dm['vz']

    # masses
    m_dm = dm['mass']*mass_fac
    Ndm = len(dm)
    print(f"Number of dark matter particles: {Ndm:.3g}")

    # particle IDs
    id_dm = np.arange(1, Ndm+1)

    
    ## Disk particles (old stars) (particle type 2)

    # position
    x_d = disk['x']
    y_d = disk['y']
    z_d = disk['z']

    # velocities
    vx_d = disk['vx']
    vy_d = disk['vy']
    vz_d = disk['vz']

    # masses
    m_d = disk['mass']*mass_fac
    Ndisk = len(disk)
    print(f"Number of disk particles: {Ndisk:.3g}")

    # particle IDs
    id_d = np.arange(1, Ndisk+1)

    
    ## Bulge particles (old stars) (particle type 3)

    # position
    x_b = bulge['x']
    y_b = bulge['y']
    z_b = bulge['z']

    # velocities
    vx_b = bulge['vx']
    vy_b = bulge['vy']
    vz_b = bulge['vz']

    # masses
    m_b = bulge['mass']*mass_fac
    Nbulge = len(bulge)
    print(f"Number of bulge particles: {Nbulge:.3g}")

    # particle IDs
    id_b = np.arange(1, Nbulge+1)


    ## open HDF5 file for writing
    file = h5py.File(filename,'w') 

    # NOTE: The following code is copied from the GIZMO public repository, mostly including Phil's comments:

    # set particle number header field [very important!]
    npart = np.array([Ngas,Ndm,Ndisk,Nbulge,0,0])

    # now we make the Header - the formatting here is peculiar, for historical (GADGET-compatibility) reasons
    h = file.create_group("Header");
    # here we set all the basic numbers that go into the header
    # (most of these will be written over anyways if it's an IC file; the only thing we actually *need* to be 'correct' is "npart")
    h.attrs['NumPart_ThisFile'] = npart; # npart set as above - this in general should be the same as NumPart_Total, it only differs 
                                         #  if we make a multi-part IC file. with this simple script, we aren't equipped to do that.
    h.attrs['NumPart_Total'] = npart; # npart set as above
    h.attrs['NumPart_Total_HighWord'] = 0*npart; # this will be set automatically in-code (for GIZMO, at least)
    h.attrs['MassTable'] = np.zeros(6); # these can be set if all particles will have constant masses for the entire run. however since 
                                        # we set masses explicitly by-particle this should be zero. that is more flexible anyways, as it 
                                        # allows for physics which can change particle masses 

    ## all of the parameters below will be overwritten by whatever is set in the run-time parameterfile if
    ##   this file is read in as an IC file, so their values are irrelevant. they are only important if you treat this as a snapshot
    ##   for restarting. Which you shouldn't - it requires many more fields be set. But we still need to set some values for the code to read
 
    h.attrs['Time'] = 0.0;  # initial time
    h.attrs['Redshift'] = 0.0; # initial redshift
    h.attrs['BoxSize'] = 1.0; # box size
    h.attrs['NumFilesPerSnapshot'] = 1; # number of files for multi-part snapshots
    h.attrs['Omega0'] = 1.0; # z=0 Omega_matter
    h.attrs['OmegaLambda'] = 0.0; # z=0 Omega_Lambda
    h.attrs['HubbleParam'] = 1.0; # z=0 hubble parameter (small 'h'=H/100 km/s/Mpc)
    h.attrs['Flag_Sfr'] = 0; # flag indicating whether star formation is on or off
    h.attrs['Flag_Cooling'] = 0; # flag indicating whether cooling is on or off
    h.attrs['Flag_StellarAge'] = 0; # flag indicating whether stellar ages are to be saved
    h.attrs['Flag_Metals'] = 0; # flag indicating whether metallicity are to be saved
    h.attrs['Flag_Feedback'] = 0; # flag indicating whether some parts of springel-hernquist model are active
    h.attrs['Flag_DoublePrecision'] = 0; # flag indicating whether ICs are in single/double precision
    h.attrs['Flag_IC_Info'] = 0; # flag indicating extra options for ICs

    ## ok, that ends the block of 'useless' parameters
    
    # Now, the actual data!
    #   These blocks should all be written in the order of their particle type (0,1,2,3,4,5)
    #   If there are no particles of a given type, nothing is needed (no block at all)
    #   PartType0 is 'special' as gas. All other PartTypes take the same, more limited set of information in their ICs
    
    # start with particle type zero. first (assuming we have any gas particles) create the group 
    p = file.create_group("PartType0")
    # now combine the xyz positions into a matrix with the correct format
    q=np.zeros((Ngas,3))
    q[:,0]=x_g
    q[:,1]=y_g
    q[:,2]=z_g
    # write it to the 'Coordinates' block
    p.create_dataset("Coordinates",data=q)
    # similarly, combine the xyz velocities into a matrix with the correct format
    q=np.zeros((Ngas,3))
    q[:,0]=vx_g
    q[:,1]=vy_g
    q[:,2]=vz_g
    # write it to the 'Velocities' block
    p.create_dataset("Velocities",data=q)
    # write particle ids to the ParticleIDs block
    p.create_dataset("ParticleIDs",data=id_g)
    # write particle masses to the Masses block
    p.create_dataset("Masses",data=m_g)
    # write internal energies to the InternalEnergy block
    p.create_dataset("InternalEnergy",data=u_g)
    # combine the xyz magnetic fields into a matrix with the correct format
    q=np.zeros((Ngas,3))
    q[:,0]=bx_g
    q[:,1]=by_g
    q[:,2]=bz_g
    # write magnetic fields to the MagneticField block. note that this is unnecessary if the code is compiled with 
    #   MAGNETIC off. however, it is not a problem to have the field there, even if MAGNETIC is off, so you can 
    #   always include it with some dummy values and then use the IC for either case
    p.create_dataset("MagneticField",data=q)

    # PartType1 for this IC (dark matter)
    p = file.create_group("PartType1")
    q=np.zeros((Ndm,3))
    q[:,0]=x_dm
    q[:,1]=y_dm
    q[:,2]=z_dm
    p.create_dataset("Coordinates",data=q)
    q=np.zeros((Ndm,3))
    q[:,0]=vx_dm
    q[:,1]=vy_dm
    q[:,2]=vz_dm
    p.create_dataset("Velocities",data=q)
    p.create_dataset("ParticleIDs",data=id_dm)
    p.create_dataset("Masses",data=m_dm)

    # PartType2 for this IC (disk particles)
    p = file.create_group("PartType2")
    q=np.zeros((Ndisk,3))
    q[:,0]=x_d
    q[:,1]=y_d
    q[:,2]=z_d
    p.create_dataset("Coordinates",data=q)
    q=np.zeros((Ndisk,3))
    q[:,0]=vx_d
    q[:,1]=vy_d
    q[:,2]=vz_d
    p.create_dataset("Velocities",data=q)
    p.create_dataset("ParticleIDs",data=id_d)
    p.create_dataset("Masses",data=m_d)

    # PartType3 for this IC (bulge particles)
    p = file.create_group("PartType3")
    q=np.zeros((Nbulge,3))
    q[:,0]=x_b
    q[:,1]=y_b
    q[:,2]=z_b
    p.create_dataset("Coordinates",data=q)
    q=np.zeros((Nbulge,3))
    q[:,0]=vx_b
    q[:,1]=vy_b
    q[:,2]=vz_b
    p.create_dataset("Velocities",data=q)
    p.create_dataset("ParticleIDs",data=id_b)
    p.create_dataset("Masses",data=m_b)

    # no PartType4 for this IC
    # no PartType5 for this IC

    file.close()


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', help='output hdf5 filename')
    feature_parser = parser.add_mutually_exclusive_group(required=False)
    feature_parser.add_argument('--magnetic-fields', dest='magnetic_fields', action='store_true')
    feature_parser.add_argument('--no-magnetic-fields', dest='magnetic_fields', action='store_false')
    parser.set_defaults(magnetic_fields=True)
    args = parser.parse_args()

    print(f"Magnetic fields: {args.magnetic_fields}")
    make_agora_IC(args.filename, magnetic_fields=args.magnetic_fields)
