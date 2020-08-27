import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.cm import colors

if __name__ == '__main__':

    ## plot cooling curve for WNM/CNM in typical conditions

    pc_in_cm = 3.0e18 # cm == 1 pc

    gamma_eos = 5./3.
    boltzmann_cgs = 1.380658e-16 # erg K^{-1}
    m_H = 1.6733e-24    # g
    mu = 0.62111795     # dimensionless mean molecular weight
    mean_molecular_weight = mu * m_H  # fully ionized solar metallicity gas (not quite right!)
    c_v = (boltzmann_cgs / mean_molecular_weight) * (1.0 / (gamma_eos - 1.0))

    def coolrate_compton(T, n_H, X_e=1e-3):
        fac = 5.65e-36
        T_cmb = 2.2755
        n_e = X_e * n_H
        fac1 = (n_e / n_H) * (T - T_cmb) / n_H
        rate = fac * fac1
        return rate

    def coolrate_dust(T, Z=1.0, T_dust=30.):
        fac = 1.12e-32
        fac1 = (T-T_dust) * (T)**(1/2)
        fac2 = (1 - 0.8*np.exp(-75/T)) * Z
        rate = fac * fac1 * fac2
        return rate

    def coolrate_photoionization(T, n_H, X_H=0.75):
        n_He = n_H * (1 - X_H)/X_H # compute number density of He in cm^{-3}

        sigma_H = 6.0e-18
        l = 4.4 * (T/1e4)**(-0.173) * pc_in_cm
        tau = sigma_H * n_H * l
        f = (1 + 0.1) * np.exp(-tau)

        log10_eps_H = -24.6 - np.exp(50.*(-1.05))
        log10_eps_He = log10_eps_H - 0.0366
        eps_H = 10**log10_eps_H
        eps_He = 10**log10_eps_He
        term1 = eps_H
        term2 = eps_He * (n_He / n_H)
        rate = -f * (term1 + term2) / n_H
        return rate

    def coolrate_cold(T, n_H, Z=1.0):
        # T is in Kelvins
        # n_H is the number density of hydrogen in cm^{-3}
        # Z is the average metallicity in units of solar metallicity

        sigma_H = 6.0e-18
        l = 4.4 * (T/1e4)**(-0.173) * pc_in_cm
        tau = sigma_H * n_H * l

        fac = 2.896e-26
        fac1 = 1.0 / ( (T/125.215)**(-4.9202) + (T/1349.86)**(-1.7288) + (T/6450.06)**(-0.3075) )
        fac2 = (1 + Z)/(1 + 0.00143*n_H) * (1 - np.exp(-tau))
        fac3 = (0.001 + (0.1*n_H)/(1+n_H) + (0.09*n_H)/(1 + 0.1*n_H) + (Z**2)/(1+n_H))
        fac4 = np.exp(-(T/1.58e5)**2)
        rate = fac * fac1 * fac2 * fac3 * fac4
        return rate

    def coolrate_cosmicrays(n_H, Z=1.0, X_H=0.75, X_e=1e-3):
        n_e = X_e*n_H
        f_CR = n_H / (0.01 + n_H)
        e_CR = 9.0e-12 * f_CR
        fac = -1.0e-16
        fac1 = (0.98 + 1.65 * (n_e/n_H) * X_H) * e_CR / n_H
        rate = fac * fac1
        return rate

    def coolrate_photoelectric(T, n_H, Z=1.0, X_e=1e-3):
        n_e = X_e*n_H
        e_photon = 1.0  # Habing FUV intensity
        x = (e_photon * T**(0.5)) / (0.5 * (n_e/n_H) * n_H)
        fac = -1.3e-24
        fac1 = (e_photon / n_H) * Z
        fac2 = 0.049/(1 + (x/1925)**(0.73)) + (0.037*(T/1e4)**(0.7))/(1+(x/5e3))
        rate = fac * fac1 * fac2
        return rate

    lognHmin = -3.0
    lognHmax = 4.0
    logTmin = 1.0
    logTmax = 4.0
    logT_range = np.linspace(logTmin, logTmax, 100)
    lognH_range = np.linspace(lognHmin, lognHmax, 100)
    T_range = 10.**(logT_range)
    nH_range = 10.**(lognH_range)
    nH_grid, T_grid = np.meshgrid(nH_range, T_range)

    coolrates  = coolrate_cold(T_grid, nH_grid)
    coolrates += coolrate_cosmicrays(nH_grid)
    coolrates += coolrate_photoelectric(T_grid, nH_grid)
    #coolrates += coolrate_photoionization(T_grid, nH_grid)
    coolrates += coolrate_dust(T_grid)
    coolrates += coolrate_compton(T_grid, nH_grid)

    def plot_coolrates(nH_grid, T_grid, coolrates, filename="coolrate.png"):
        coolrates_per_vol = nH_grid**2 * coolrates
        rho = m_H * nH_grid
        E_thermal = rho * c_v * T_grid
        cooltime = np.abs( E_thermal / coolrates_per_vol / 3.15e7 )

        ## plot cooltime
        fig, ax = plt.subplots(figsize=(6,4))
        x = lognH_range
        y = logT_range
        p = ax.imshow(cooltime,
                        cmap='viridis',
                        extent=[x.min(), x.max(), y.min(), y.max()],
                        interpolation='nearest',
                        origin='lower',
                        aspect='auto',
                        norm=colors.LogNorm(vmin=1e2, vmax=1e9))

        fig.colorbar(p)
        ax.set_title("cooling time [yr]")
        ax.set_xlabel(r"\log n_H (H cm$^{-3}$)")
        ax.set_ylabel(r"\log T (Kelvin)")
        fig.tight_layout()
        fig_dpi = 300
        plt.savefig(filename, dpi=fig_dpi)
        plt.close()

    plot_coolrates(nH_grid, T_grid, coolrates, filename="coolrates_neutralgas_noUVB.png")

    coolrates += coolrate_photoionization(T_grid, nH_grid)
    plot_coolrates(nH_grid, T_grid, coolrates, filename="coolrates_neutralgas.png")