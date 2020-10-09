from Meshoid import Meshoid
from load_from_snapshot import load_from_snapshot
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import colors

# physical constants
boltzmann_cgs = 1.380658e-16 # erg K^{-1}
m_H = 1.6733e-24    # g
gamma = 5./3.   # assumed constant

# conversion factor from specific energy in code units to specific energy in cgs
unittime_cgs = 3.08568e16    # s  (0.976 Gyr in s)
unitlength_cgs = 3.085678e21 # cm (1 kpc in cm)
unitmass_cgs = 1.989e43      # g  (1e10 Msun in g)
unitvelocity_cgs = unitlength_cgs / unittime_cgs
unitenergypermass_cgs = unitlength_cgs**(2) * unittime_cgs**(-2) # erg g^{-1}
unitdensity_cgs = unitmass_cgs * unitlength_cgs**(-3) # g cm^{-3}
unitdensity_per_H = unitdensity_cgs / m_H
unitenergydensity_cgs = unitmass_cgs * unitvelocity_cgs**(2) * unitlength_cgs**(-3)
unitbfield_cgs = np.sqrt(4*np.pi*unitenergydensity_cgs)

print(f"unit B field in gauss = {unitbfield_cgs}")

# plot resolution
fig_dpi = 300

def load_hydro_data(filename):
    pdata = {}
    for field in "Masses", "Coordinates", "SmoothingLength", "Velocities", "InternalEnergy", "ElectronAbundance", "MagneticField":
        pdata[field] = load_from_snapshot(field, 0, filename)

    pdata["Metallicity"] = load_from_snapshot("Metallicity", 0, filename)
    return pdata

def load_stars(filename):
    star_coords = load_from_snapshot("Coordinates", 4, filename)
    if star_coords is None:
        return None
    else:
        return star_coords

def compute_temperature(pdata):
    # ONLY if COOL_LINE_METALS is enabled
    #metal_mass_frac = pdata["Metallicity"][:,0] # total metallicity (mass fraction)
    #He_mass_frac = pdata["Metallicity"][:,1]
    #y_Helium = He_mass_frac
    #z_metals = metal_mass_frac

    # otherwise, assume solar metallicity
    y_Helium = 0.23
    z_metals = 0.02

    x_H = 1.0 - y_Helium - z_metals

    H_term = x_H
    He_term = y_Helium/4.0
    z_term = z_metals/2.0   # this term is ignored

    e_term_fullyionized = x_H + y_Helium*(2)/4.0
    mu_fullyionized = 1.0 / (H_term + He_term + e_term_fullyionized)
    #print(f"fully-ionized mean molecular weight (dimensionless): {mu_fullyionized}")
    #print(f"maximium electron abundance: {np.max(pdata['ElectronAbundance'])}")

    InternalEnergy = pdata["InternalEnergy"] * unitenergypermass_cgs
    if pdata['ElectronAbundance'] is None:
        return None

    e_term = pdata["ElectronAbundance"]

    mu = 1.0 / (H_term + He_term + e_term)
    mean_molecular_weight = mu * m_H
    T = (mean_molecular_weight / boltzmann_cgs) * (gamma-1) * InternalEnergy
    return T

def apply_radius_cut(pdata, T, rmax=40.):
    coords = pdata["Coordinates"]
    radius_cut = np.sum(coords*coords,axis=1) < (rmax**2)

    ndata = {}
    ndata['Coordinates']     = coords[radius_cut]
    ndata['Masses']          = pdata['Masses'][radius_cut]
    ndata['SmoothingLength'] = pdata['SmoothingLength'][radius_cut]
    ndata['Velocities']      = pdata['Velocities'][radius_cut]
    ndata['MagneticField']   = pdata['MagneticField'][radius_cut]
    if T is not None:
        temp = T[radius_cut]
    else:
        temp = None
        
    return ndata, temp

def save_phase_plot(input_dens, temp, filename):
    ## make phase plot! (temperature vs n_H)
    dens = input_dens * unitdensity_per_H

    lognHmin = -7.5
    lognHmax = 5.5
    logTmin = np.log10(3.0)
    logTmax = 8.5
    
    h = plt.hist2d(np.log10(dens), np.log10(temp),
                    bins=100, norm=colors.LogNorm(),
                    range = [[lognHmin,lognHmax], [logTmin,logTmax]])

    plt.colorbar(label=r'proportional to mass') # unclear what units of this are
    plt.xlabel(r'$\log_{10}$ density ($n_{H}$ [cm$^{-3}$])')
    plt.ylabel(r'$\log_{10}$ temperature (K)')
    plt.savefig(filename, dpi=fig_dpi)
    plt.close()

def plot_stars_on_axis(ax, star_coords, rmax=10., s=0.05):
    if star_coords is not None:
        x = star_coords[:,0]
        y = star_coords[:,1]
        r = np.sqrt(x**2 + y**2)
        ax.scatter(x[r < rmax], y[r < rmax], s=s, color='black')

def save_slice_plot(mesh, field, filename, colorbar_label="",
                    bfield=None,
                    star_coords=None, rmax=10., plane='z', vmin=10., vmax=1.0e7):
    res = 800
    x = y = np.linspace(-rmax,rmax,res)
    X, Y = np.meshgrid(x, y)

    #sigma_gas_msun_pc2 = 1e4 * M.SurfaceDensity(M.m,center=np.array([0,0,0]),size=40.,res=res)
    #density_slice_nHcgs =300 * M.Slice(M.Density(),center=np.array([0,0,0]),size=40.,res=res)
    slice = mesh.Slice(field,center=np.array([0,0,0]),size=2*rmax,res=res,plane=plane)

    fig,ax = plt.subplots(figsize=(6,6))
    p = ax.imshow(slice, cmap='viridis',
                    extent=[x.min(), x.max(), y.min(), y.max()],
                    interpolation='nearest',
                    origin='lower',
                    aspect='equal',
                    norm=colors.LogNorm(vmin=vmin, vmax=vmax))

    bfield_res = 20
    bfield_slice = mesh.Slice(bfield,
                              center=np.array([0,0,0]),
                              size=2*rmax,
                              res=bfield_res)
    
    x = y = np.linspace(-rmax, rmax, bfield_res)
    bfield_x = bfield_slice[:,:,0]
    bfield_y = bfield_slice[:,:,1]
    bfield_z = bfield_slice[:,:,2]
    bfield_norm = np.sqrt( bfield_x**2 + bfield_y**2 + bfield_z**2 )
    bfield_x /= bfield_norm
    bfield_y /= bfield_norm
    # plot B-fields
    ax.quiver(x, y, bfield_x, bfield_y, angles='xy', pivot='middle')

    plot_stars_on_axis(ax, star_coords, rmax=rmax)

    ax.set_aspect('equal')
    #fig.colorbar(p,label=r"$\Sigma_{gas}$ $(\rm M_\odot\,pc^{-2})$")
    #fig.colorbar(p, label=r"$n_{H}$ (cm$^{-3}$)")
    fig.colorbar(p, label=colorbar_label)
    ax.set_xlabel("X (kpc)")
    ax.set_ylabel("Y (kpc)")
    plt.savefig(filename, dpi=fig_dpi)
    plt.close()

def save_density_projection_plot(mesh, filename, star_coords=None,
                                bfield=None,
                                rmax=10.,
                                colorbar_label=r"$\Sigma_{gas}$ $(\rm M_\odot\,pc^{-2})$"):
    res = 1000
    x = y = np.linspace(-rmax,rmax,res)
    X, Y = np.meshgrid(x, y)

    sigma_gas_msun_pc2 = 1e4 * mesh.SurfaceDensity(mesh.m,
                                    center=np.array([0,0,0]),size=2*rmax,res=res)

    sigma_gas_msun_pc2 = np.swapaxes(sigma_gas_msun_pc2, 0, 1)

    fig,ax = plt.subplots(figsize=(6,6))
    p = ax.imshow(sigma_gas_msun_pc2, cmap='viridis',
                    extent=[x.min(), x.max(), y.min(), y.max()],
                    interpolation='nearest',
                    origin='lower',
                    aspect='equal',
                    norm=colors.LogNorm(vmin=1.0, vmax=1.0e3))

    bfield_res = 20
    bfield_slice = mesh.Slice(bfield,
                              center=np.array([0,0,0]),
                              size=2*rmax,
                              res=bfield_res)
    
    x = y = np.linspace(-rmax, rmax, bfield_res)
    bfield_x = bfield_slice[:,:,0]
    bfield_y = bfield_slice[:,:,1]
    bfield_z = bfield_slice[:,:,2]
    bfield_norm = np.sqrt( bfield_x**2 + bfield_y**2 + bfield_z**2 )
    bfield_x /= bfield_norm
    bfield_y /= bfield_norm
    # plot B-fields
    ax.quiver(x, y, bfield_x, bfield_y, angles='xy', pivot='middle')

    plot_stars_on_axis(ax, star_coords, rmax=rmax)

    ax.set_aspect('equal')
    fig.colorbar(p, label=colorbar_label)
    ax.set_xlabel("X (kpc)")
    ax.set_ylabel("Y (kpc)")
    plt.savefig(filename, dpi=fig_dpi)
    plt.close()

def compute_mesh(pdata):
    return Meshoid(pdata["Coordinates"], pdata["Masses"], pdata["SmoothingLength"])
