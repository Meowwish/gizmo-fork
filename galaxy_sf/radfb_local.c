#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "../allvars.h"
#include "../proto.h"
#include "../kernel.h"

/* independent re-implementation of photoionization feedback by Armillotta et al. */

#ifdef GALSF_PHOTOIONIZATION
static inline void downheap2 (double* data1, double* data2, double* data3, int* data4, int* data5, const size_t N, size_t k) {
    
    double v1 = data1[k];
    double v2 = data2[k];
    double v3 = data3[k];
    int v4 = data4[k];
    int v5 = data5[k];

    while (k<=N/2) {
        size_t j = 2 * k;
        if (j < N && data1[j] < data1[(j + 1)]) j++;
        if (!(v1 < data1[j]))  break;
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


void sort (double *data1, double* data2, double* data3, int* data4, int* data5, const size_t n) {
    
    if (n==0) return;
    size_t N = n - 1;
    size_t k = N / 2;
    k++;     
    do {
        k--;
        downheap2(data1, data2, data3, data4, data5, N, k);
    }
    while (k > 0);

    while (N > 0) {
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
    double *Distance = (double *)malloc(N_gas * sizeof(double));
    double *IonRate = (double *)malloc(N_gas * sizeof(double));
    double *Tini = (double *)malloc(N_gas * sizeof(double));
    int *ParticleNum = (int *)malloc(N_gas * sizeof(int));
    int *Tag_HIIregion = (int *)malloc(N_gas * sizeof(int));

    double Tfin = 1e4;                                       // if gas is photoionized, assume it is heated to this temperature
    double molw_i = 4.0 / (8 - 5 * (1 - HYDROGEN_MASSFRAC)); /* mean molecular weight assuming full ionization (approximately ~0.6) */
    // case B recombination coefficient (approximate)
    double beta = 3e-13; //cm**3 s*-1

    // Wake stars after one star time step //
    for (int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if (P[i].Type != 4)
            continue;
        if (P[i].Mass <= 0)
            continue;
        for (int j = 0; j < N_gas; j++)
            if (SphP[j].photo_star == P[i].ID && SphP[j].HIIregion == 1)
            {
                SphP[j].HIIregion = 0;
                SphP[j].photo_subtime = 0;
            }
    }

    for (int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        if (P[i].Type != 4)
            continue;
        if (P[i].Mass <= 0)
            continue;

        /* assume a fidicual number of 10^49 ionizing photons per 100 solar masses */
        double N_photons_per_100Msun = 1.0e49; // s^-1

        double star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
        if(star_age >= 0.005) {
            N_photons_per_100Msun = 0.0; // assume no ionizing photons after 5 Myr
        }

        double stellar_mass_Msun = P[i].Mass * (All.UnitMass_in_g / SOLAR_MASS);
        double N_photons = N_photons_per_100Msun * (stellar_mass_Msun / 100.);
        if (N_photons <= 0)
            continue;

#ifdef GALSF_PHOTOIONIZATION_DEBUGGING
        double n_H = 10.; // lower bound for HII region mean density
        double r1_approx = pow(3.0 * N_photons / (4.0 * M_PI * n_H * n_H * beta), 1. / 3.);
        printf("[Photoionization] Q [photons/sec]: %g\n", N_photons);
        printf("\tApproximate upper bound on size of HII region [code units]: %g\n", r1_approx / All.UnitLength_in_cm);
#endif

        double gas_timestep;
        double star_timestep = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;

        // FIXME: this *only* loops over the local particles!!
        for (int j = 0; j < N_gas; j++) /* loop over the gas block */
        {
            ParticleNum[j] = j;
            Tag_HIIregion[j] = SphP[j].HIIregion;
            double distx = P[i].Pos[0] - P[j].Pos[0];
            double disty = P[i].Pos[1] - P[j].Pos[1];
            double distz = P[i].Pos[2] - P[j].Pos[2];
            Distance[j] = sqrt(distx * distx + disty * disty + distz * distz);

            double Rhob = SphP[j].Density * UNIT_DENSITY_IN_CGS;
            double Mb = P[j].Mass * All.UnitMass_in_g;

            // compute temperature of fluid element
            Tini[j] = CallGrackle(SphP[j].InternalEnergy, SphP[j].Density, 0, SphP[j].Ne, j, 2);

            // dimensionless mean molecular weight of fluid element (prior to photoionization)
            // N.B. SphP[j].InternalEnergy is the *specific* (per unit mass) internal energy
            double molw_n = Tini[j] * BOLTZMANN / (EOS_GAMMA - 1) / (SphP[j].InternalEnergy * UNIT_ENERGY_IN_CGS / All.UnitMass_in_g) / PROTONMASS;

            IonRate[j] = HYDROGEN_MASSFRAC * beta * Rhob * Mb / (2 * PROTONMASS * PROTONMASS * molw_n * molw_i);
        }

        //Cycle to sort particles by increasing distance
        sort(Distance, IonRate, Tini, ParticleNum, Tag_HIIregion, N_gas);

        // FIXME: this *only* loops over local particles!!
        int j;
        for (j = 0; j < N_gas; j++) /* loop over the gas block */
        {
            if (Tag_HIIregion[j] == 1)
                continue; // The particle belongs to another HII region
            if (Tini[j] >= Tfin)
                continue; //Particle already ionized
            if (IonRate[j] <= N_photons)
            {
                SphP[ParticleNum[j]].InternalEnergy = BOLTZMANN * Tfin / ((EOS_GAMMA - 1) * molw_i * PROTONMASS) * All.UnitMass_in_g / UNIT_ENERGY_IN_CGS;
                SphP[ParticleNum[j]].InternalEnergyPred = BOLTZMANN * Tfin / ((EOS_GAMMA - 1) * molw_i * PROTONMASS) * All.UnitMass_in_g / UNIT_ENERGY_IN_CGS;
                SphP[ParticleNum[j]].HIIregion = 1;
                gas_timestep = (P[ParticleNum[j]].TimeBin ? (1 << P[ParticleNum[j]].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
                SphP[ParticleNum[j]].photo_subtime = round(star_timestep / gas_timestep);
                SphP[ParticleNum[j]].photo_star = P[i].ID;
                N_photons -= IonRate[j];
            }
            else
            {
                double Prandom = get_random_number(ThisTask);
                if (IonRate[j] / N_photons > Prandom)
                {
                    SphP[ParticleNum[j]].InternalEnergy = BOLTZMANN * Tfin / ((EOS_GAMMA - 1) * molw_i * PROTONMASS) * All.UnitMass_in_g / UNIT_ENERGY_IN_CGS;
                    SphP[ParticleNum[j]].InternalEnergyPred = BOLTZMANN * Tfin / ((EOS_GAMMA - 1) * molw_i * PROTONMASS) * All.UnitMass_in_g / UNIT_ENERGY_IN_CGS;
                    SphP[ParticleNum[j]].HIIregion = 1;
                    gas_timestep = (P[ParticleNum[j]].TimeBin ? (1 << P[ParticleNum[j]].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a;
                    SphP[ParticleNum[j]].photo_subtime = round(star_timestep / gas_timestep);
                    SphP[ParticleNum[j]].photo_star = P[i].ID;
                    N_photons -= IonRate[j];
                }
            }
            //P[i].Feedback_timestep = DMIN(P[i].Feedback_timestep,fbtime);
            if (N_photons <= 0)
                break;
        }

#ifdef GALSF_PHOTOIONIZATION_DEBUGGING
        int jmax = (j < N_gas) ? j : (N_gas - 1);
        double r1 = Distance[jmax];
        printf("\tActual size of HII region [code units]: %g\n", r1);
#endif
    }

    free(IonRate);
    free(Distance);
    free(Tini);
    free(ParticleNum);
    free(Tag_HIIregion);
}

#endif // GALSF_PHOTOIONIZATION