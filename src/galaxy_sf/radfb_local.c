#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

#ifdef SLUG
#include "slug_feedback.hpp"
#endif

/* independent re-implementation of photoionization feedback by Armillotta et al. */

#ifdef GALSF_PHOTOIONIZATION
static inline void downheap2(double *data1, double *data2, double *data3, int *data4, int *data5, const size_t N, size_t k)
{

    double v1 = data1[k];
    double v2 = data2[k];
    double v3 = data3[k];
    int v4 = data4[k];
    int v5 = data5[k];

    while (k <= N / 2)
    {
        size_t j = 2 * k;
        if (j < N && data1[j] < data1[(j + 1)])
            j++;
        if (!(v1 < data1[j]))
            break;
        data1[k] = data1[j];
        data2[k] = data2[j];
        data3[k] = data3[j];
        data4[k] = data4[j];
        data5[k] = data5[j];
        k = j;
    }
    data1[k] = v1;
    data2[k] = v2;
    data3[k] = v3;
    data4[k] = v4;
    data5[k] = v5;
}

void sort(double *data1, double *data2, double *data3, int *data4, int *data5, const size_t n)
{

    if (n == 0)
        return;
    size_t N = n - 1;
    size_t k = N / 2;
    k++;
    do
    {
        k--;
        downheap2(data1, data2, data3, data4, data5, N, k);
    } while (k > 0);

    while (N > 0)
    {
        /* first swap the elements */
        double tmp;

        tmp = data1[0];
        data1[0] = data1[N];
        data1[N] = tmp;

        tmp = data2[0];
        data2[0] = data2[N];
        data2[N] = tmp;

        tmp = data3[0];
        data3[0] = data3[N];
        data3[N] = tmp;

        int tmp2;

        tmp2 = data4[0];
        data4[0] = data4[N];
        data4[N] = tmp2;

        tmp2 = data5[0];
        data5[0] = data5[N];
        data5[N] = tmp2;

        N--;
        downheap2(data1, data2, data3, data4, data5, N, 0);
    }
}

void compute_photoionization(void)
{
    double *Distance = (double *) mymalloc("Distance", N_gas * sizeof(double));
    double *IonRate = (double *) mymalloc("IonRate", N_gas * sizeof(double));
    double *Tini = (double *) mymalloc("Tini", N_gas * sizeof(double));
    int *ParticleNum = (int *) mymalloc("ParticleNum", N_gas * sizeof(int));
    int *Tag_HIIregion = (int *) mymalloc("Tag_HIIregion", N_gas * sizeof(int));

    // if gas is photoionized, assume it is heated to this temperature
    const double Tfin = 1.0e4;
    // mean molecular weight assuming full ionization (approximately ~0.6)
    const double molw_i = 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC));
    // case B recombination coefficient (approximate)
    const double beta = 3.0e-13; // cm**3 s*-1

#if 0
    // Wake gas particles after one star time step
    for (int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if (P[i].Type != 4)
        {
            continue;
        }
        if (P[i].Mass <= 0)
        {
            continue;
        }
        for (int j = 0; j < N_gas; j++) {
            if (SphP[j].photo_star == P[i].ID && SphP[j].HIIregion == 1) // race condition (photo_star could have moved to a different processor)
            {
                SphP[j].HIIregion = 0;
                SphP[j].photo_subtime = 0;
            }
        }
    }
#endif

    // reset HII region tags on *all* local gas particles
    for (int j = 0; j < N_gas; j++)
    {
        SphP[j].HIIregion = 0;
    }

    // loop over *all* local star particles
    //for (int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    for (int i = 0; i < NumPart; i++)
    {
        if (P[i].Type != 4)
        {
            continue;
        }
        if (P[i].Mass <= 0)
        {
            continue;
        }

#ifdef SLUG
        double N_photons = slugComputeIonizingPhotons(i); // compute number of ionizing photons via SLUG
#else
        // *without* SLUG: assume a fidicual number of 10^49 ionizing photons per 100 solar masses
        const double N_photons_per_100Msun = 1.0e49; // s^-1
        const double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);

        // assume no ionizing photons after 5 Myr
        if (star_age >= 0.005)
        {
            N_photons_per_100Msun = 0.0;
        }

        const double stellar_mass_Msun = P[i].Mass * UNIT_MASS_IN_SOLAR;
        double N_photons = N_photons_per_100Msun * (stellar_mass_Msun / 100.);
#endif // SLUG

        if (N_photons <= 0.)
        {
            continue;
        }

#ifdef GALSF_PHOTOIONIZATION_DEBUGGING
        const double n_H = 100.;                                                                   // approximate lower bound for HII region mean density
        const double r1_approx = pow(3.0 * N_photons / (4.0 * M_PI * n_H * n_H * beta), 1. / 3.); // cm
        const double cm_in_parsec = 3.085678e18;
        printf("[Photoionization] Q [photons/sec/(100 Msun)] = %g\n", N_photons / (P[i].Mass * UNIT_MASS_IN_SOLAR / 100.));
        printf("\tApproximate upper bound on size of HII region = %g pc; ", r1_approx / cm_in_parsec);
#endif

        //const double star_timestep = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;

        // loop over all local gas particles
        for (int j = 0; j < N_gas; j++) /* loop over the gas block */
        {
            ParticleNum[j] = j;
            Tag_HIIregion[j] = SphP[j].HIIregion;
            const double dx = P[i].Pos[0] - P[j].Pos[0];
            const double dy = P[i].Pos[1] - P[j].Pos[1];
            const double dz = P[i].Pos[2] - P[j].Pos[2];
            Distance[j] = sqrt(dx * dx + dy * dy + dz * dz);

            const double Rhob = SphP[j].Density * UNIT_DENSITY_IN_CGS;
            const double Mb = P[j].Mass * UNIT_MASS_IN_CGS;

            // compute temperature of fluid element
            Tini[j] = CallGrackle(SphP[j].InternalEnergy, SphP[j].Density, 0, SphP[j].Ne, j, 2);

            // dimensionless mean molecular weight of fluid element (prior to photoionization)
            // N.B. SphP[j].InternalEnergy is the internal energy *per unit mass*
            const double molw_n = Tini[j] * BOLTZMANN / (EOS_GAMMA - 1) / (SphP[j].InternalEnergy * UNIT_ENERGY_IN_CGS / UNIT_MASS_IN_CGS) / PROTONMASS;

            IonRate[j] = HYDROGEN_MASSFRAC * beta * Rhob * Mb / (2 * PROTONMASS * PROTONMASS * molw_n * molw_i);
        }

        // sort particles by increasing distance
        sort(Distance, IonRate, Tini, ParticleNum, Tag_HIIregion, N_gas);

        // loop over all local gas particles
        int jmax = (N_gas - 1);
        for (int j = 0; j < N_gas; j++) /* loop over the gas block */
        {
            if (Tag_HIIregion[j] == 1)
            {
                continue; // The particle belongs to another HII region
            }

            if (Tini[j] >= Tfin)
            {
                continue; // Particle already ionized
            }

            const double e_ion = BOLTZMANN * Tfin / ((EOS_GAMMA - 1) * molw_i * PROTONMASS) * UNIT_MASS_IN_CGS / UNIT_ENERGY_IN_CGS;
            //const double gas_timestep = (P[ParticleNum[j]].TimeBin ? (1 << P[ParticleNum[j]].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;

            if (IonRate[j] <= N_photons)
            {
                SphP[ParticleNum[j]].InternalEnergy = e_ion;
                SphP[ParticleNum[j]].InternalEnergyPred = e_ion;
                SphP[ParticleNum[j]].HIIregion = 1;

                //SphP[ParticleNum[j]].photo_subtime = round(star_timestep / gas_timestep);
                //SphP[ParticleNum[j]].photo_star = P[i].ID;

                N_photons -= IonRate[j];
            }
            else
            {
                const double Prandom = get_random_number(ThisTask);

                if (IonRate[j] / N_photons > Prandom)
                {
                    SphP[ParticleNum[j]].InternalEnergy = e_ion;
                    SphP[ParticleNum[j]].InternalEnergyPred = e_ion;
                    SphP[ParticleNum[j]].HIIregion = 1;

                    //SphP[ParticleNum[j]].photo_subtime = round(star_timestep / gas_timestep);
                    //SphP[ParticleNum[j]].photo_star = P[i].ID;

                    N_photons -= IonRate[j];
                }
            }

            if (N_photons <= 0)
            {
                jmax = j;
                break;
            }
        }

#ifdef GALSF_PHOTOIONIZATION_DEBUGGING
        const double r1 = Distance[jmax];
        printf("actual size of HII region = %g pc.\n", r1 * UNIT_LENGTH_IN_PC);
#endif
    }

    // free temporary arrays
    // NOTE: this *MUST* be done in exactly the reverse order that they are allocated above!
    myfree(Tag_HIIregion);
    myfree(ParticleNum);
    myfree(Tini);
    myfree(IonRate);
    myfree(Distance);
}

#endif // GALSF_PHOTOIONIZATION
