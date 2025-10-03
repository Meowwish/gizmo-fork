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

/* Routines for mechanical feedback/enrichment models: stellar winds, supernovae, etc
 * This file was written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 */

#if defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)

// Track the number of SNIa that have occurred in the simulation
static double last_SNIa_output_time = 0.0;
static int total_SNIa_since_last_output = 0;
static int total_SNIa_all_time = 0;

void determine_where_SNe_occur(void)
{
    if (All.Time <= 0) return;

    double dt, star_age, npossible, nhosttotal, ntotal, ptotal, dtmean, rmean;
    npossible = nhosttotal = ntotal = ptotal = dtmean = rmean = 0;
    double mpi_npossible, mpi_nhosttotal, mpi_ntotal, mpi_ptotal, mpi_dtmean, mpi_rmean;
    mpi_npossible = mpi_nhosttotal = mpi_ntotal = mpi_ptotal = mpi_dtmean = mpi_rmean = 0;
    double total_old_stellar_mass = 0.0;
    int SNIa_this_step = 0;
    
    // calculate the total mass of old stellar particles
    double local_old_stellar_mass = 0.0;
    for(int i = 0; i < NumPart; i++)
    {
        if((P[i].Type == 4 && P[i].StellarAge < 0 && P[i].Mass > 0) || (P[i].Type == 2))
        {
            local_old_stellar_mass += P[i].Mass * UNIT_MASS_IN_SOLAR;
        }
    }
    
    MPI_Reduce(&local_old_stellar_mass, &total_old_stellar_mass, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Bcast(&total_old_stellar_mass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // 添加保护，防止total_old_stellar_mass为0
    if(total_old_stellar_mass <= 0.0)
    {
        total_old_stellar_mass = 1.0; // 使用默认值1.0以避免除以零错误
    }

    //printf("P2DEBUG Total mass of old stellar particles: %g Msun\n", total_old_stellar_mass);
    
    // loop over particles //
    const double sn_loop_begin_walltime = MPI_Wtime();
    for (int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        P[i].SNe_ThisTimeStep=0;

#if defined(SINGLE_STAR_SINK_DYNAMICS)
        if (P[i].Type == 0)
        {
            continue;
        } // any non-gas type is eligible to be a 'star' here
#else
        if (All.ComovingIntegrationOn)
        {
            if (P[i].Type != 4)
            {
                continue;
            }
        } // in cosmological simulations, 'stars' have particle type=4
        if (All.ComovingIntegrationOn == 0)
        {
            if ((P[i].Type < 2) || (P[i].Type > 4))
            {
                continue;
            }
        }  // in non-cosmological sims, types 2,3,4 are valid 'stars'
#endif
        if(P[i].Mass<=0) {continue;}
        dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
#ifdef BH_INTERACT_ON_GAS_TIMESTEP
	if(P[i].Type == 5){dt = P[i].dt_since_last_gas_search;}
#endif
        if(dt<=0) {continue;} // no time, no events
        star_age = evaluate_stellar_age_Gyr(P[i].StellarAge);
        if (star_age <= 0)
        {
            continue;
        } // unphysical age, no events

        // StellarAge<0，old stellar particles (not formed from SLUG) According to allvars.h, StellarAge equals to StellarFormationTime
        if ((P[i].Type == 4 && P[i].StellarAge < 0) || (P[i].Type == 2))// Meow_OldStarCheck
        {
            // Calculate SNIa events for old stellar particles
            double SNIa_rate = 5000.0; // SNIa events rate - 5000 Myr^-1
            double NSM_rate= 500.0; // NSM events rate - 500 Myr^-1
            double particle_mass_solar = P[i].Mass * UNIT_MASS_IN_SOLAR; // Stellar mass in solar mass
            double dt_Myr = dt * UNIT_TIME_IN_MYR; // Time step in Myr
            
            // Expected SNIa number in this time step = event_rate * time_step * mass_fraction
            double expected_SNIa = SNIa_rate * dt_Myr * (particle_mass_solar / total_old_stellar_mass);
            double expected_NSM = NSM_rate * dt_Myr * (particle_mass_solar / total_old_stellar_mass);
            
            // Calculate the number of SNIa events for this particle
            int num_SNIa = 0;
            int num_NSM = 0;
            if (expected_SNIa > 1.0) {
                // Expected SNIa number is large, take the integer part
                num_SNIa = (int)(expected_SNIa + get_random_number(i));
            } else if (expected_SNIa > 0.0) {
                // Expected SNIa number is small, use a random number to decide
                if (get_random_number(i) < expected_SNIa) num_SNIa = 1;
            }

            if (expected_NSM > 1.0) {
                // Expected NSM number is large, take the integer part
                num_NSM = (int)(expected_NSM + get_random_number(i));
            } else if (expected_NSM > 0.0) {
                // Expected NSM number is small, use a random number to decide
                if (get_random_number(i) < expected_NSM) num_NSM = 1;
            }
            
            // Update the particle properties
            if (num_SNIa > 0) {
                // DEBUG: Print SNIa event information
                //printf("SNIa_DEBUG: Particle ID=%d generating %d SNIa events. Mass=%g Msun, dt=%g Myr, expected=%g\n", 
                //       P[i].ID, num_SNIa, particle_mass_solar, dt_Myr, expected_SNIa);
                //fflush(stdout);
                
                P[i].SNe_ThisTimeStep = num_SNIa;
                ntotal += num_SNIa;
                nhosttotal++;
                
                // Count the SNIa events for this step
                SNIa_this_step += num_SNIa;
            }            
            
            // Update the particle properties
            if (num_NSM > 0) {
                // DEBUG: Print SNIa event information
                //printf("SNIa_DEBUG: Particle ID=%d generating %d SNIa events. Mass=%g Msun, dt=%g Myr, expected=%g\n", 
                //       P[i].ID, num_SNIa, particle_mass_solar, dt_Myr, expected_SNIa);
                //fflush(stdout);
                
                P[i].SNe_ThisTimeStep = num_NSM * 10;
            }
            
            npossible++;
            continue;
        }
        npossible++; // it is possible for a SN event to occur

//if (P[i].SNe_ThisTimeStep > 0) {
//    printf("SNe_ThisTimeStep_BeforeFB: %.1f, Cumulative: %.1f\n", P[i].SNe_ThisTimeStep, P[i].SNe_Cumulative);
//}
#ifdef SLUG
        slugComputeSNFeedback(i);
// Meow: check if the SLUG feedback is working
//if (P[i].SNe_ThisTimeStep > 0) {
//    printf("SNe_ThisTimeStep_AfterFB: %.1f, Cumulative: %.1f\n", P[i].SNe_ThisTimeStep, P[i].SNe_Cumulative);
//}

#else  // *without* SLUG
        double RSNe = mechanical_fb_calculate_eventrates(i, dt);
        rmean += RSNe;
        ptotal += RSNe * (P[i].Mass * UNIT_MASS_IN_SOLAR) * (dt * UNIT_TIME_IN_MYR);
#endif // SLUG

#ifdef GALSF_SFR_IMF_SAMPLING
        if (P[i].IMF_NumMassiveStars > 0)
        {
            P[i].IMF_NumMassiveStars = DMAX(0, P[i].IMF_NumMassiveStars - P[i].SNe_ThisTimeStep);
        } // lose an O-star for every SNe //
#endif

        if (P[i].SNe_ThisTimeStep > 0)
        {
            ntotal += P[i].SNe_ThisTimeStep;
            nhosttotal++;
        }
        dtmean += dt;
    } // for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i]) //

#ifdef SLUG_DEBUG_SN_RATE
    // collect SNe data
    int mpi_SNIa_this_step = 0;
    MPI_Reduce(&SNIa_this_step, &mpi_SNIa_this_step, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (ThisTask == 0) {
        total_SNIa_since_last_output += mpi_SNIa_this_step;
        total_SNIa_all_time += mpi_SNIa_this_step;
        
        // check if we need to output the statistics
        double current_time_Myr = All.Time * UNIT_TIME_IN_MYR;
        double last_output_time_Myr = last_SNIa_output_time * UNIT_TIME_IN_MYR;
        
        // if the current time is greater than the last output time, or if this is the first output
        if (floor(current_time_Myr * 10) > floor(last_output_time_Myr * 10) || last_SNIa_output_time == 0) {
            // output the statistics
            printf("SNIa_STATS: Time=%g Myr, SNIa events in last %.2f Myr: %d, Total events: %d, Rate: %.2f per Myr\n", 
                   current_time_Myr, 
                   current_time_Myr - last_output_time_Myr,
                   total_SNIa_since_last_output,
                   total_SNIa_all_time,
                   (current_time_Myr > last_output_time_Myr) ? 
                       total_SNIa_since_last_output / (current_time_Myr - last_output_time_Myr) : 0);
            fflush(stdout);
            
            // reset the last output time and count
            last_SNIa_output_time = All.Time;
            total_SNIa_since_last_output = 0;
        }
    }
    
    // broadcast the last SNIa output time and total SNIa count to all tasks
    MPI_Bcast(&last_SNIa_output_time, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&total_SNIa_all_time, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // end timing measurement
    //const double sn_loop_end_walltime = MPI_Wtime();
    // compute total mass in star particles (P[i].Type == 4)
    double thisTaskStellarMass = 0.; // solar masses
    double thisTaskNumberSNe = 0.; // dimensionless

    for(int i=0; i < NumPart; ++i)
    {
      if (P[i].Type != 4) {
	    continue;
      }

      if (P[i].Mass <= 0.) {
	    continue;
      }

      const double mass = P[i].Mass * UNIT_MASS_IN_SOLAR;
      const double numberSNe = P[i].SNe_Cumulative;

      thisTaskStellarMass += mass;
      thisTaskNumberSNe += numberSNe;
    }
    
    double totalStellarMass; // solar masses
    double totalNumberSNe; // dimensionless
    MPI_Reduce(&thisTaskStellarMass, &totalStellarMass, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&thisTaskNumberSNe, &totalNumberSNe, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (ThisTask == 0)
    {
      const double SNe_per_100Msun = totalNumberSNe / (totalStellarMass / 100.0);
      std::cout << "[SLUG] Total stellar mass = " << totalStellarMass << " Msun\n";
      std::cout << "[SLUG] Total number of SNe = " << totalNumberSNe << "\n";
      std::cout << "[SLUG] SNe per 100 Msun = " << SNe_per_100Msun << std::endl;
    }
#endif // SLUG_DEBUG_SN_RATE

    // end timing measurement
    const double sn_loop_end_walltime = MPI_Wtime();

#ifdef SLUG_DEBUG_PERFORMANCE
    double mpi_snloop_begin_time;
    double mpi_snloop_end_time;

    MPI_Reduce(&sn_loop_begin_walltime, &mpi_snloop_begin_time, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    MPI_Reduce(&sn_loop_end_walltime, &mpi_snloop_end_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if (ThisTask == 0)
    {
        const double slug_elapsed_time = mpi_snloop_end_time - mpi_snloop_begin_time;
        slug_total_elapsed_time += slug_elapsed_time;

        std::cout << "[SLUG] Processed SLUG objects in "
                  << slug_elapsed_time
                  << " seconds.\n"
                  << "[SLUG] SLUG accounts for "
                  << 100. * (slug_total_elapsed_time / CPUThisRun)
                  << "% of overall runtime.\n";
    }
#endif // SLUG_DEBUG_PERFORMANCE

    MPI_Reduce(&dtmean, &mpi_dtmean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&rmean, &mpi_rmean, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ptotal, &mpi_ptotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&nhosttotal, &mpi_nhosttotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&ntotal, &mpi_ntotal, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&npossible, &mpi_npossible, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(ThisTask == 0)
    {
        if (mpi_ntotal > 0 && mpi_nhosttotal > 0 && mpi_dtmean > 0)
            if (mpi_npossible > 0)
            {
                mpi_dtmean /= mpi_npossible;
                mpi_rmean /= mpi_npossible;
                fprintf(FdSneIIHeating, "%lg %g %g %g %g %g %g \n", All.Time, mpi_npossible, mpi_nhosttotal, mpi_ntotal, mpi_ptotal, mpi_dtmean, mpi_rmean);
            }
        if (All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
        {
            mpi_dtmean /= mpi_npossible; mpi_rmean /= mpi_npossible;
            fprintf(FdSneIIHeating, "%lg %g %g %g %g %g %g \n", All.Time,mpi_npossible,mpi_nhosttotal,mpi_ntotal,mpi_ptotal,mpi_dtmean,mpi_rmean); fflush(FdSneIIHeating);
        }
    } // if(ThisTask == 0) //

} // void determine_where_SNe_occur() //

#endif // defined(GALSF_FB_MECHANICAL) || defined(GALSF_FB_THERMAL)


// The rest of this file is for mechanical feedback only

#ifdef GALSF_FB_MECHANICAL


int addFB_evaluate_active_check(int i, int fb_loop_iteration);
int addFB_evaluate_active_check(int i, int fb_loop_iteration)
{
    if(P[i].Type <= 1) { return 0; }
    if(PPP[i].Hsml <= 0) { return 0; }
    if(PPP[i].NumNgb <= 0) { return 0; }
#ifdef BH_INTERACT_ON_GAS_TIMESTEP
    if(P[i].Type == 5 && !P[i].do_gas_search_this_timestep) { return 0; }
#endif

    // SNII, SNIa, NSM events always activate the particle
    if(P[i].SNe_ThisTimeStep > 0) {
        if(fb_loop_iteration < 0 || fb_loop_iteration == 0) { return 1; }
    }
    // SLUG AGB activate the particle
    else if(P[i].Type == 4 && P[i].StellarAge > 30 / UNIT_TIME_IN_MYR && TimeBinActive[P[i].TimeBin]) { // not dynamically active in this time bin, don't activate for chemical feedback!!!
        if(fb_loop_iteration < 0 || fb_loop_iteration == 0) { return 1; }
    }
    // old AGB activate the particle
    else if(P[i].Type == 4 && P[i].StellarAge < 0) {
        if(!TimeBinActive[P[i].TimeBin]) { return 0; } // not dynamically active in this time bin, don't activate for chemical feedback!!!
        if(fb_loop_iteration <= 0) { 
            // use the particle ID and current simulation time to create a pseudo-random selection
            unsigned int hash = (P[i].ID + (unsigned int)(All.Time * 1000)) * 2654435761u;
            if(hash % 20 == 0) {  // 5% of particles are selected
                return 1;
            }
        }
    }
    // stellar disc activate the particle
    else if(P[i].Type == 2) {
        if(!TimeBinActive[P[i].TimeBin]) { return 0; } // not dynamically active in this time bin, don't activate for chemical feedback!!!
        if(fb_loop_iteration <= 0) { 
            // use the particle ID and current simulation time to create a pseudo-random selection
            unsigned int hash = (P[i].ID + (unsigned int)(All.Time * 1000)) * 2654435761u;
            if(hash % 20 == 0) {  // 5% of particles are selected
                return 1;
            }
         }
    }
    return 0; // not active
}
/*
int addFB_evaluate_active_check(int i, int fb_loop_iteration)
{
    if(P[i].Type <= 1) {return 0;} // note quantities used here must -not- change in the loop [hence not using mass here], b/c can change offsets for return from different processors, giving a negative mass and undefined behaviors
    if(PPP[i].Hsml <= 0) {return 0;}
    if(PPP[i].NumNgb <= 0) {return 0;}
#ifdef BH_INTERACT_ON_GAS_TIMESTEP
    if(P[i].Type == 5 && !P[i].do_gas_search_this_timestep) {return 0;}
#endif
    if(P[i].SNe_ThisTimeStep>0) {if(fb_loop_iteration<0 || fb_loop_iteration==0) {return 1;}}
    
}
*/
#define CORE_FUNCTION_NAME addFB_evaluate /* name of the 'core' function doing the actual inter-neighbor operations. this MUST be defined somewhere as "int CORE_FUNCTION_NAME(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)" */
#define INPUTFUNCTION_NAME particle2in_addFB    /* name of the function which loads the element data needed (for e.g. broadcast to other processors, neighbor search) */
#define OUTPUTFUNCTION_NAME out2particle_addFB  /* name of the function which takes the data returned from other processors and combines it back to the original elements */
#define CONDITIONFUNCTION_FOR_EVALUATION if(addFB_evaluate_active_check(i,loop_iteration)) /* function for which elements will be 'active' and allowed to undergo operations. can be a function call, e.g. 'density_is_active(i)', or a direct function call like 'if(P[i].Mass>0)' */
#include "../system/code_block_xchange_initialize.h" /* pre-define all the ALL_CAPS variables we will use below, so their naming conventions are consistent and they compile together, as well as defining some of the function calls needed */

// define kernel structure (purely for convenience, will hold variables below) //
struct kernel_addFB {double dp[3], r, wk, dwk, hinv, hinv3, hinv4;};

struct OUTPUT_STRUCT_NAME
{
  MyFloat M_coupled;
  MyFloat Area_weighted_sum[AREA_WEIGHTED_SUM_ELEMENTS];
#ifdef SN_INJECTED_MOMENTUM_ACCOUNTING // Momentum/energy summation for each event
  MyFloat injected_radial_momentum;
  MyFloat injected_thermal_energy;
  MyFloat pnorm;
  MyFloat pvec[3]; // face vectors
#endif // SN_INJECTED_MOMENTUM_ACCOUNTING
}
*DATARESULT_NAME, *DATAOUT_NAME;


void particle2in_addFB(struct addFB_evaluate_data_in_ *in, int i, int loop_iteration)
{
    // pre-assign various values that will be used regardless of feedback physics //
    int k; for(k=0;k<3;k++) {in->Pos[k] = P[i].Pos[k]; in->Vel[k] = P[i].Vel[k];}
    double heff = PPP[i].Hsml / PPP[i].NumNgb; in->V_i = heff*heff*heff; in->Hsml = PPP[i].Hsml;
#ifdef METALS
    for(k=0;k<NUM_METAL_SPECIES;k++) {in->yields[k]=0.0;}
#endif
    for(k=0;k<AREA_WEIGHTED_SUM_ELEMENTS;k++) {in->Area_weighted_sum[k] = P[i].Area_weighted_sum[k];}
    in->Msne = 0; in->unit_mom_SNe = 0; in->SNe_v_ejecta = 0;
    if((P[i].DensAroundStar <= 0)||(P[i].Mass <= 0)) {return;} // events not possible [catch for mass->0]
    if(loop_iteration < 0) {in->Msne=P[i].Mass; in->unit_mom_SNe=1.e-4; in->SNe_v_ejecta=1.0e-4; return;} // weighting loop

    // Check the SN event is SNII or SNIa
    if((P[i].Type == 4 && P[i].StellarAge < 0 && P[i].SNe_ThisTimeStep > 0 && P[i].SNe_ThisTimeStep < 10 && loop_iteration == 0)||(P[i].Type == 2 && P[i].SNe_ThisTimeStep > 0 && P[i].SNe_ThisTimeStep < 10 && loop_iteration == 0))
    {
        // MODIFIED: Don't call particle2in_addFB_fromstars for old stars to avoid SLUG errors
        // Instead set SNIa parameters directly
        in->Msne = 1.4 / UNIT_MASS_IN_SOLAR; // Typical SNIa ejecta mass (1 Msun)
        in->SNe_v_ejecta = sqrt(2.*(1.e51/UNIT_ENERGY_IN_CGS)/in->Msne); // Typical SNIa ejecta velocity (10,000 km/s)
#ifdef METALS
        // SNIa yields in solar masses (Seitenzahl et al. 2012, Table 2, Z=0.02 column) TEST FINISHED, NO NEEDED ANYMORE
        double snia_yields[31] = {
            0.0,      // metallicity[0] for cooling
            0.0,      // SNII   C12     1   SLUG
            0.0,      // SNII   N14     2   SLUG
            0.0,      // SNII   O16     3   SLUG
            0.0,      // SNII   Mg24    4   SLUG
            0.0,      // SNII   S32     5   SLUG
            0.0,      // WR     C12     6   SLUG
            0.0,      // WR     N14     7   SLUG
            0.0,      // WR     O16     8   SLUG
            0.0,      // WR     Mg24    9   SLUG
            0.0,      // WR     S32     10  SLUG
            0.0,      // AGB    C12     11  SLUG
            0.0,      // AGB    N14     12  SLUG
            0.0,      // AGB    O16     13  SLUG
            0.0,      // AGB    Ba138   14  SLUG
            0.0,      // AGB    Ce140   15  SLUG
            6.22E-01, // SNIa   Fe56    16
            6.90E-02, // SNIa   Ni58    17
            3.21E-06, // SNIa   N14     18
            2.84E-01, // SNIa   Si28    19
            1.11E-01, // SNIa   S32     20
            0.0,      // old_AGB        21
            0.0,      // old_AGB        22
            0.0,      // old_AGB        23
            0.0,      // old_AGB        24
            0.0,      // old_AGB        25
            0.0,      // NSM            26
            0.0,      // NSM            27
            0.0,      // NSM            28
            0.0,      // NSM            29
            0.0       // NSM            30
        };
        // Set metal yields (multiply by number of SNIa events and scale by mass)
        //for(k=16; k < 21; k++) {
        //    in->yields[k] += snia_yields[k] * P[i].SNe_ThisTimeStep;
        //}
        in->yields[40] += 1.4 * P[i].SNe_ThisTimeStep;
#endif
        in->unit_mom_SNe = in->Msne * in->SNe_v_ejecta;
        return;
    }

    if((P[i].Type == 4 && P[i].StellarAge < 0 && P[i].SNe_ThisTimeStep > 9 && loop_iteration == 0)||(P[i].Type == 2 && P[i].SNe_ThisTimeStep > 9 && loop_iteration == 0))
    {
        // MODIFIED: Don't call particle2in_addFB_fromstars for old stars to avoid SLUG errors
        // Instead set SNIa parameters directly
        in->Msne = 0.18 / UNIT_MASS_IN_SOLAR; // Typical NSM ejecta mass (1 Msun)
        in->SNe_v_ejecta = sqrt(2.*(1.e51/UNIT_ENERGY_IN_CGS)/in->Msne); // Typical NSM ejecta velocity
        in->yields[41] += 0.18 * P[i].SNe_ThisTimeStep/10;
        in->unit_mom_SNe = in->Msne * in->SNe_v_ejecta;
        return;
    }
    
    // Add the old AGB feedback
    if((P[i].Type == 4 && P[i].StellarAge < 0 && loop_iteration == 0 && P[i].SNe_ThisTimeStep == 0)||(P[i].Type == 2 && loop_iteration == 0 && P[i].SNe_ThisTimeStep == 0))
    {
#ifdef METALS
        double particle_mass_solar = P[i].Mass * UNIT_MASS_IN_SOLAR;
        double dt = GET_PARTICLE_TIMESTEP_IN_PHYSICAL(i);
        double dt_Myr = dt * UNIT_TIME_IN_MYR; // Time step in Myr
        //printf("AGB total mass loss: %g, particle mass: %g", particle_mass_solar * dt_Myr * 3E-06 / UNIT_MASS_IN_SOLAR, particle_mass_solar);
        in->Msne += particle_mass_solar * dt_Myr * 1.93E-05 / UNIT_MASS_IN_SOLAR * 20;
        //for(k=21; k < 26; k++) {
        //    in->yields[k] += agb_yields[k] * dt_Myr * particle_mass_solar;
            //if (k == 21) {printf("\nAGB yields: C12 = %g\n", in->yields[k]);}
        //}
        in->yields[42] += particle_mass_solar * dt_Myr * 1.93E-05 * 20;
        return;
#endif
    }

    particle2in_addFB_fromstars(in,i,loop_iteration); // subroutine that actually deals with the assignment of feedback properties
    in->unit_mom_SNe = in->Msne * in->SNe_v_ejecta;
}

void out2particle_addFB(struct OUTPUT_STRUCT_NAME *out, int i, int mode, int loop_iteration)
{
    if(P[i].Mass > 0)
    {
        if(loop_iteration < 0)
        {
            int k=0, kmin=0, kmax=7; if(loop_iteration == -1) {kmin=kmax; kmax=AREA_WEIGHTED_SUM_ELEMENTS;}
#ifdef GALSF_USE_SNE_ONELOOP_SCHEME
            kmin=0; kmax=AREA_WEIGHTED_SUM_ELEMENTS;
#endif
            for(k=kmin;k<kmax;k++) {ASSIGN_ADD(P[i].Area_weighted_sum[k], out->Area_weighted_sum[k], mode);}
        } else {
            P[i].Mass -= out->M_coupled;
            if ((P[i].Mass < 0) || (isnan(P[i].Mass)))
            {
                P[i].Mass = 0;
            }
#ifdef SN_INJECTED_MOMENTUM_ACCOUNTING // Momentum/energy accounting
            const double momentum_cgs = out->injected_radial_momentum * UNIT_MASS_IN_CGS * UNIT_VEL_IN_CGS;
            const double momentum_per_Msun_cgs = (momentum_cgs / SOLAR_MASS);   // cm/s
            const double momentum_per_Msun_kms = momentum_per_Msun_cgs / 1.0e5; // km/s
            const double energy_cgs = out->injected_thermal_energy * UNIT_ENERGY_IN_CGS;

            P[i].SNe_InjectedMomentumThisStep += momentum_per_Msun_kms;
            P[i].SNe_InjectedThermalEnergyThisStep += energy_cgs;
            P[i].SNe_pnorm += out->pnorm;
            for(int k=0; k<3; ++k) {
                P[i].SNe_pvec[k] += out->pvec[k];
            }
#endif // SN_INJECTED_MOMENTUM_ACCOUNTING
        }
    }
}


/* here we have the subroutine that is the work center of this module. does the key calculations over neighbors and actually couples the relevant feedback quantities */
/*!   -- this subroutine [in all versions] writes important conservative variables to shared memory [updating the neighbor values]:
        we need to protect these writes for openmp (thread safety) below. the read-in values are themselves modified, so we need to protect -both- the read and write operations in all cases */

#ifdef GALSF_USE_SNE_ONELOOP_SCHEME

#ifdef SLUG
#error "Do not use SNE_ONELOOP_SCHEME with SLUG. Unphysical things may happen!"
#endif // SLUG

int addFB_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int startnode, numngb_inbox, listindex = 0, j, k, n;
    double u,r2,h2,kernel_zero,wk,dM,dP,E_coupled,dP_sum,dP_boost_sum;
    struct kernel_addFB kernel; struct addFB_evaluate_data_in_ local; struct OUTPUT_STRUCT_NAME out;
    memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME));
    kernel_main(0.0,1.0,1.0,&kernel_zero,&wk,-1); wk=0;
    if(mode == 0) {particle2in_addFB(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];} /* Load the data for the particle injecting feedback */
    if(local.Msne<=0) {return 0;} // no SNe for the origin particle! nothing to do here //
    if(local.Hsml<=0) {return 0;} // zero-extent kernel, no particles //
    h2 = local.Hsml*local.Hsml; kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);
    double unitlength_in_kpc=UNIT_LENGTH_IN_KPC * All.cf_atime, density_to_n=All.cf_a3inv*UNIT_DENSITY_IN_NHCGS, unit_egy_SNe = 1.0e51/UNIT_ENERGY_IN_CGS; // some units (just used below, but handy to define for clarity) //


    // now define quantities that will be used below //
    double Esne51; Esne51 = 0.5*local.SNe_v_ejecta*local.SNe_v_ejecta*local.Msne / unit_egy_SNe;
    double RsneKPC, RsneKPC_0; RsneKPC=0.; RsneKPC_0=(0.0284/unitlength_in_kpc) * pow(1+Esne51,0.286); //Cioffi: weak external pressure
    double r2max_phys = 2.0/unitlength_in_kpc; r2max_phys *= r2max_phys; // no super-long-range effects allowed! (of course this is arbitrary in code units) //

    /* Now start the actual FB computation for this particle */
    if(mode == 0) {startnode = All.MaxPart;} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;} /* root node & node opening */
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_pairs_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb_inbox < 0) {return -2;}

            E_coupled = dP_sum = dP_boost_sum = 0;
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if(P[j].Type != 0) {continue;} // require a gas particle //
                
                double Mass_j, InternalEnergy_j, rho_j, Vel_j[3]; // initialize holders for the local variables that might change below
                #pragma omp atomic read
                Mass_j = P[j].Mass; // this can get modified below, so we need to read it thread-safe now
                
                // quick block of checks to make sure it's actually worth continuing!
                if(Mass_j <= 0) continue; // require the particle has mass //
                for(k=0; k<3; k++) {kernel.dp[k] = local.Pos[k] - P[j].Pos[k];}
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2], 1); // find the closest image in the given box size  //
                r2=0; for(k=0;k<3;k++) {r2 += kernel.dp[k]*kernel.dp[k];}
                if(r2<=0) {continue;} // same particle //
                double h2j = PPP[j].Hsml * PPP[j].Hsml;
                if((r2>h2)&&(r2>h2j)) {continue;} // outside kernel (in both 'directions') //
                if(r2 > r2max_phys) {continue;} // outside long-range cutoff //
                kernel.r = sqrt(r2); if(kernel.r <= 0) {continue;}
                
                // calculate kernel quantities //
                #pragma omp atomic read
                rho_j = SphP[j].Density;
                u = kernel.r * kernel.hinv;
                double hinv_j = 1./PPP[j].Hsml, hinv3_j = hinv_j*hinv_j*hinv_j; /* note these lines and many below assume 3D sims! */
                double wk_j = 0, dwk_j = 0, u_j = kernel.r * hinv_j, hinv4_j = hinv_j*hinv3_j, V_j = Mass_j / rho_j;
                if(u<1) {kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 1);} else {kernel.dwk=kernel.wk=0;}
                if(u_j<1) {kernel_main(u_j, hinv3_j, hinv4_j, &wk_j, &dwk_j, 1);} else {wk_j=dwk_j=0;}
                if(local.V_i<0 || isnan(local.V_i)) {local.V_i=0;}
                if(V_j<0 || isnan(V_j)) {V_j=0;}
                double sph_area = fabs(local.V_i*local.V_i*kernel.dwk + V_j*V_j*dwk_j); // effective face area //
                wk = 0.5 * (1 - 1/sqrt(1 + sph_area / (M_PI*kernel.r*kernel.r))); // corresponding geometric weight //
                if((wk <= 0)||(isnan(wk))) continue; // no point in going further, there's no physical weight here
                double wk_vec[AREA_WEIGHTED_SUM_ELEMENTS]={0}, wk_tmp=0;
                wk_vec[0] = wk;
                wk_tmp=wk*kernel.dp[0]/kernel.r; if(kernel.dp[0]>0) {wk_vec[1]=wk_tmp;} else {wk_vec[2]=wk_tmp;}
                wk_tmp=wk*kernel.dp[1]/kernel.r; if(kernel.dp[1]>0) {wk_vec[3]=wk_tmp;} else {wk_vec[4]=wk_tmp;}
                wk_tmp=wk*kernel.dp[2]/kernel.r; if(kernel.dp[2]>0) {wk_vec[5]=wk_tmp;} else {wk_vec[6]=wk_tmp;}

                // if loop_iteration==-1, this is a pre-calc loop to get the relevant weights for coupling //
                if(loop_iteration < 0)
                {
                    for(k=0;k<AREA_WEIGHTED_SUM_ELEMENTS;k++) out.Area_weighted_sum[k] += wk_vec[k];
                    continue;
                }
                // NOW do the actual feedback calculation //
                double wk_norm = 1. / (MIN_REAL_NUMBER + fabs(local.Area_weighted_sum[0])); // normalization for scalar weight sum
                wk *= wk_norm; // this way wk matches the value summed above for the weighting //
                if((wk <= 0)||(isnan(wk))) continue;

                // ok worth initializing other variables we will use below
                #pragma omp atomic read
                InternalEnergy_j = SphP[j].InternalEnergy; // this can get modified below, so we need to read it thread-safe now
                for(k=0;k<3;k++) {
                    #pragma omp atomic read
                    Vel_j[k] = P[j].Vel[k]; // this can get modified below, so we need to read it thread-safe now
                }
                double InternalEnergy_j_0 = InternalEnergy_j, Mass_j_0 = Mass_j, rho_j_0 = rho_j, Vel_j_0[3]; for(k=0;k<3;k++) {Vel_j_0[k]=Vel_j[k];} // save initial values to use below
#ifdef METALS
                double Metallicity_j[NUM_METAL_SPECIES], Metallicity_j_0[NUM_METAL_SPECIES];
                for(k=0;k<NUM_METAL_SPECIES;k++) {
                    #pragma omp atomic read
                    Metallicity_j[k] = P[j].Metallicity[k]; // this can get modified below, so we need to read it thread-safe now
                    Metallicity_j_0[k] = Metallicity_j[k]; // save initial values to  use below
                }
#endif
                
                /* define initial mass and ejecta velocity in this 'cone' */
                double v_bw[3]={0}, e_shock=0, pnorm = 0, pvec[3]={0};
                for(k=0; k<3; k++)
                {
                    double q; q = 0; int i1=2*k+1, i2=i1+1;
                    double q_i1 = fabs(local.Area_weighted_sum[i1]);
                    double q_i2 = fabs(local.Area_weighted_sum[i2]);
                    if((q_i1>MIN_REAL_NUMBER)&&(q_i2>MIN_REAL_NUMBER))
                    {
                        double rr = q_i2/q_i1;
                        double rr2 = rr * rr;
                        if(wk_vec[i1] != 0)
                        {
                            q += wk_norm * wk_vec[i1] * sqrt(0.5*(1.0+rr2));
                        } else {
                            q += wk_norm * wk_vec[i2] * sqrt(0.5*(1.0+1.0/rr2));
                        }
                    } else {
                        q += wk_norm * (wk_vec[i1] + wk_vec[i2]);
                    }
                    pvec[k] = -q;
                    pnorm += pvec[k]*pvec[k];
                }
                pnorm = sqrt(pnorm);

                wk = pnorm; // this (vector norm) is the new 'weight function' for our purposes
                dM = wk * local.Msne;

                /* now, add contribution from relative star-gas particle motion to shock energy */
                for(k=0;k<3;k++)
                {
                    v_bw[k] = local.SNe_v_ejecta*pvec[k]/pnorm + (local.Vel[k]-Vel_j[k])/All.cf_atime;
                    e_shock += v_bw[k]*v_bw[k];
                }
                double mj_preshock, dM_ejecta_in, massratio_ejecta, mu_j;
                mj_preshock = Mass_j;
                dM_ejecta_in = dM;
                massratio_ejecta = dM_ejecta_in / (dM_ejecta_in + Mass_j);
                mu_j = Mass_j / (dM + Mass_j);
                e_shock *= pnorm * 0.5*local.Msne * mu_j;

                if((wk <= 0)||(isnan(wk))) continue;

                RsneKPC = RsneKPC_0;
                double n0 = rho_j*density_to_n;
                /* this is tedious, but is a fast approximation (essentially a lookup table) for the -0.429 power above */
                if(n0 < 1.e-3) {RsneKPC *= 19.4;} else {
                    if(n0 < 1.e-2) {RsneKPC *= 1.9 + 23./(1.+333.*n0);} else {
                        if(n0 < 1.e-1) {RsneKPC *= 0.7 + 8.4/(1.+33.3*n0);} else {
                            if(n0 < 1) {RsneKPC *= 0.08 + 3.1/(1.+2.5*n0);} else {
                                if(n0 < 10) {RsneKPC *= 0.1 + 1.14/(1.+0.333*n0);} else {
                                    if(n0 < 100) {RsneKPC *= 0.035 + 0.43/(1.+0.0333*n0);} else {
                                        if(n0 < 1000) {RsneKPC *= 0.017 + 0.154/(1.+0.00333*n0);} else {
                                            if(n0 < 1.e4) {RsneKPC *= 0.006 + 0.057/(1.+0.000333*n0);} else {
                                                RsneKPC *= pow(n0, -0.429); }}}}}}}}


                /* below expression is again just as good a fit to the simulations, and much faster to evaluate */
                double z0 = Metallicity_j[0]/All.SolarAbundances[0];
                if(z0 < 0.01) {RsneKPC *= 2.0;} else {
                    if(z0 < 1) {RsneKPC *= 0.93 + 0.0615 / (0.05 + 0.8*z0);} else {RsneKPC *= 0.8 + 0.4 / (1 + z0);}}
                /* calculates cooling radius given density and metallicity in this annulus into which the ejecta propagate */

                /* if coupling radius > R_cooling, account for thermal energy loss in the post-shock medium:
                 from Thornton et al. thermal energy scales as R^(-6.5) for R>R_cool */
                double r_eff_ij = sqrt(r2) - Get_Particle_Size(j);
                if(r_eff_ij > RsneKPC) {e_shock *= RsneKPC*RsneKPC*RsneKPC/(r_eff_ij*r_eff_ij*r_eff_ij);}

                /* now we have the proper energy to couple */
                E_coupled += e_shock;

                /* inject actual mass from mass return */
                int couple_anything_but_scalar_mass_and_metals = 1; // key to indicate whether or not we actually need to do the next set of steps beyond pure scalar mass+metal couplings //
                if(P[j].Hsml<=0) {if(rho_j>0){rho_j*=(1+dM_ejecta_in/Mass_j);} else {rho_j=dM_ejecta_in*kernel.hinv3;}} else {rho_j+=kernel_zero*dM_ejecta_in*hinv3_j;}
                rho_j *= 1 + dM_ejecta_in/Mass_j; // inject mass at constant particle volume //
                Mass_j += dM_ejecta_in;
                out.M_coupled += dM_ejecta_in;
#if defined(METALS) /* inject metals */
                for(k=0;k<NUM_METAL_SPECIES-NUM_AGE_TRACERS;k++) {
                    double old_metal = Metallicity_j[k];
                    Metallicity_j[k] = (1 - massratio_ejecta) * Metallicity_j[k] + local.yields[k]/(local.Msne*1e10) * massratio_ejecta; // !!!!!!

                    // DEBUG: Track O16 injection
                    //if(k == 3 && local.yields[k] > 0) {
                    //    printf("SNIa_DEBUG: Metal injection to gas ID=%d: O16 from=%g to=%g, yield=%g, massratio=%g\n", 
                    //           P[j].ID, old_metal, Metallicity_j[k], local.yields[k]/(local.Msne*1e10), massratio_ejecta);
                    //    fflush(stdout);
                    //}

                }
#ifdef GALSF_FB_FIRE_AGE_TRACERS
                if(loop_iteration == 3) {for(k=NUM_METAL_SPECIES-NUM_AGE_TRACERS;k<NUM_METAL_SPECIES;k++) {Metallicity_j[k] += pnorm*local.yields[k]/Mass_j;}} // add age tracers in taking yields to mean MASS, so we can make it large without actually exchanging large masses
#ifndef GALSF_FB_FIRE_AGE_TRACERS_DISABLE_SURFACE_YIELDS
                if(loop_iteration != 3) {for(k=NUM_METAL_SPECIES-NUM_AGE_TRACERS;k<NUM_METAL_SPECIES;k++) {Metallicity_j[k]=(1-massratio_ejecta)*Metallicity_j[k] + massratio_ejecta*local.yields[k];}} // treat like any other yield when doing stellar mass exchange
#endif
#endif
#endif
                if(couple_anything_but_scalar_mass_and_metals)
                {
                    
                /* inject the post-shock energy and momentum (convert to specific units as needed first) */
                e_shock *= 1 / Mass_j;
                InternalEnergy_j += e_shock;
                /* inject momentum */
                double m_ej_input = pnorm * local.Msne;
                /* appropriate factor for the ejecta being energy-conserving inside the cooling radius (or Hsml, if thats smaller) */
                double m_cooling = 4.18879*pnorm*rho_j*RsneKPC*RsneKPC*RsneKPC;
                /* apply limiter for energy conservation */
                double mom_boost_fac = 1 + sqrt(DMIN(mj_preshock , m_cooling) / m_ej_input);
#if (defined(FLAG_NOT_IN_PUBLIC_CODE) && (FLAG_NOT_IN_PUBLIC_CODE > 2)) || defined(SINGLE_STAR_SINK_DYNAMICS) 
                if(loop_iteration > 0) {mom_boost_fac=1;} /* no unresolved PdV component for winds+r-process */
#endif
                /* save summation values for outputs */
                dP = local.unit_mom_SNe / Mass_j * pnorm;
                dP_sum += dP;
                dP_boost_sum += dP * mom_boost_fac;

                /* actually do the injection */
                double q0 = All.cf_atime * (pnorm*local.Msne/Mass_j) * mom_boost_fac;
                for(k=0; k<3; k++)
                {
                    double q = q0 * v_bw[k];
                    Vel_j[k] += q;
                }
                    
                } // couple_anything_but_scalar_mass_and_metals
                
                /* we updated variables that need to get assigned to element 'j' -- let's do it here in a thread-safe manner */
                #pragma omp atomic
                P[j].Mass += Mass_j - Mass_j_0; // finite mass update [delta difference added here, allowing for another element to update in the meantime]. done this way to ensure conservation.
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                #pragma omp atomic
                SphP[j].MassTrue += Mass_j - Mass_j_0; // finite mass update
#endif
                if(rho_j_0 > 0) {
                    #pragma omp atomic
                    SphP[j].Density *= rho_j / rho_j_0; // inject mass at constant particle volume [no need to be exactly conservative here] //
                }
                for(k=0;k<3;k++) {
                    #pragma omp atomic
                    P[j].Vel[k] += Vel_j[k] - Vel_j_0[k]; // delta-update
                    #pragma omp atomic
                    SphP[j].VelPred[k] += Vel_j[k] - Vel_j_0[k]; // delta-update
                    #pragma omp atomic
                    P[j].dp[k] += Mass_j*Vel_j[k] - Mass_j_0*Vel_j_0[k]; // discrete momentum change
                }
                #pragma omp atomic
                SphP[j].InternalEnergy += InternalEnergy_j - InternalEnergy_j_0; // delta-update
                #pragma omp atomic
                SphP[j].InternalEnergyPred += InternalEnergy_j - InternalEnergy_j_0; // delta-update
                for(k=0;k<NUM_METAL_SPECIES;k++) {
                    #pragma omp atomic
                    P[j].Metallicity[k] += Metallicity_j[k] - Metallicity_j_0[k]; // delta-update
                }
                
                
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;}}}    /* open it */
    } // while(startnode >= 0)

    /* Now collect the result at the right place */
    if(mode == 0) {out2particle_addFB(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;}
    return 0;
} // int addFB_evaluate



#else // un-protected [updated, more fixed energy-injecting SNe scheme]


int addFB_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex, int *ngblist, int loop_iteration)
{
    int startnode, numngb_inbox, listindex = 0, j, k, n;
    double u,r2,kernel_zero,wk,dM_ejecta_in,dP,E_coupled,dP_sum,dP_boost_sum;
    struct kernel_addFB kernel; struct addFB_evaluate_data_in_ local; struct OUTPUT_STRUCT_NAME out;
    memset(&out, 0, sizeof(struct OUTPUT_STRUCT_NAME));

    /* Load the data for the particle injecting feedback */
    if(mode == 0) {particle2in_addFB(&local, target, loop_iteration);} else {local = DATAGET_NAME[target];}
    if(local.Msne<=0) {return 0;} // no SNe for the origin particle! nothing to do here //
    if(local.Hsml<=0) {return 0;} // zero-extent kernel, no particles //

    // some units (just used below, but handy to define for clarity) //
    double h2 = local.Hsml*local.Hsml; kernel_main(0.0,1.0,1.0,&kernel_zero,&wk,-1); wk=0; // define the kernel zero-point value, needed to prevent some nasty behavior when no neighbors found
    kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4); // define kernel quantities
    double unitlength_in_kpc= UNIT_LENGTH_IN_KPC * All.cf_atime, density_to_n=All.cf_a3inv*UNIT_DENSITY_IN_NHCGS, unit_egy_SNe = 1.0e51/UNIT_ENERGY_IN_CGS;

    // now define quantities that will be used below //
    double psi_cool=1, psi_egycon=1, v_ejecta_eff=local.SNe_v_ejecta;
    double wk_norm = 1. / (MIN_REAL_NUMBER + fabs(local.Area_weighted_sum[0])); // normalization for scalar weight sum
    double pnorm_sum = 1./(MIN_REAL_NUMBER + fabs(local.Area_weighted_sum[10])); // re-normalization after second pass for normalized "pnorm" (should be close to ~1)
    if((local.Area_weighted_sum[0] > MIN_REAL_NUMBER) && (loop_iteration >= 0))
    {
        double vba_2_eff = pnorm_sum * local.Area_weighted_sum[7]; // phi term for energy: weighted mass-deposited KE for ejecta neighbors
        v_ejecta_eff = sqrt(local.SNe_v_ejecta*local.SNe_v_ejecta + vba_2_eff); // account for all terms to get the revised KE term here
        double beta_egycon = sqrt(pnorm_sum / local.Msne) * (1./v_ejecta_eff) * local.Area_weighted_sum[8]; // beta term for re-normalization for energy [can be positive or negative]
        double beta_cool = pnorm_sum * local.Area_weighted_sum[9]; // beta term if all particles in terminal-momentum-limit
        if(All.ComovingIntegrationOn) {if(fabs(beta_cool) < fabs(beta_egycon)) {beta_egycon = beta_cool;}}
        psi_egycon = sqrt(1. + beta_egycon*beta_egycon) - beta_egycon; // exact solution for energy equation for constant psi
        if(beta_egycon > 20.) {psi_egycon = 1./(2.*beta_egycon);} // replace with series expansion to avoid roundoff error at high beta
        if(beta_cool > 0.5) {psi_cool = 1./(2.*beta_cool);} // for cooling limit, only need upper limit to psi, all else will use less energy
    }


    double Energy_injected_codeunits = 0.5 * local.Msne * v_ejecta_eff * v_ejecta_eff;
    double Esne51 = Energy_injected_codeunits / unit_egy_SNe;
    double RsneKPC = 0., RsneKPC_3 = 0., m_cooling = 0., v_cooling = 210./UNIT_VEL_IN_KMS;
    double RsneKPC_0 = (0.0284/unitlength_in_kpc);
    int feedback_type_is_SNe = 0;
    if(loop_iteration == 0) {feedback_type_is_SNe = 1;} // assume, for now, that loop 0 represents SNe, for purposes of energy-momentum switch below //
    if(feedback_type_is_SNe == 1) // check for SNe specifically
    {
        RsneKPC_0 *= pow(1+Esne51,0.286); //SNe: using scaling from Cioffi with weak external pressure
    } else {
        RsneKPC_0 *= pow(Esne51,0.286); // ensures smooth conservation for winds and tracers as mass-loading goes to vanishingly small values
    }
    double r2max_phys = 2.0/unitlength_in_kpc; // no super-long-range effects allowed! (of course this is arbitrary in code units) //
    if(local.Hsml >= r2max_phys) {psi_egycon=DMIN(psi_egycon,1); psi_cool=DMIN(psi_cool,1);}
    r2max_phys *= r2max_phys;

    /* Now start the actual FB computation for this particle */
    if(mode == 0) {startnode = All.MaxPart;} else {startnode = DATAGET_NAME[target].NodeList[0]; startnode = Nodes[startnode].u.d.nextnode;}    /* start at root node, open it */
    while(startnode >= 0)
    {
        while(startnode >= 0)
        {
            numngb_inbox = ngb_treefind_pairs_threads(local.Pos, local.Hsml, target, &startnode, mode, exportflag, exportnodecount, exportindex, ngblist);
            if(numngb_inbox < 0) {return -2;}

            E_coupled = dP_sum = dP_boost_sum = 0;
            for(n = 0; n < numngb_inbox; n++)
            {
                j = ngblist[n]; /* since we use the -threaded- version above of ngb-finding, its super-important this is the lower-case ngblist here! */
                if(P[j].Type != 0) {continue;} // require a gas particle //
                
                double Mass_j, InternalEnergy_j, rho_j, Vel_j[3]; // initialize holders for the local variables that might change below
                #pragma omp atomic read
                Mass_j = P[j].Mass; // this can get modified below, so we need to read it thread-safe now

                // now consider a block of conditions we will use to evaluate whether its worth opening this loop at all //
                if(Mass_j <= 0) {continue;} // require the particle has mass //
                for(k=0; k<3; k++) {kernel.dp[k] = local.Pos[k] - P[j].Pos[k];}
                NEAREST_XYZ(kernel.dp[0],kernel.dp[1],kernel.dp[2],1); // find the closest image in the given box size  //
                r2=0; for(k=0;k<3;k++) {r2 += kernel.dp[k]*kernel.dp[k];}
                if(r2<=0) {continue;} // same particle //
                double h2j = PPP[j].Hsml * PPP[j].Hsml;
                if((r2>h2)&&(r2>h2j)) {continue;} // outside kernel (in both 'directions') //
                if(r2 > r2max_phys) {continue;} // outside long-range cutoff //
                kernel.r = sqrt(r2); if(kernel.r <= 0) {continue;}

                
                // calculate kernel quantities //
                #pragma omp atomic read
                rho_j = SphP[j].Density;
                u = kernel.r * kernel.hinv;
                double hinv_j = 1./PPP[j].Hsml, hinv3_j = hinv_j*hinv_j*hinv_j;
                double wk_j = 0, dwk_j = 0, u_j = kernel.r * hinv_j, hinv4_j = hinv_j*hinv3_j, V_j = Mass_j / rho_j;
                if(u<1) {kernel_main(u, kernel.hinv3, kernel.hinv4, &kernel.wk, &kernel.dwk, 1);} else {kernel.wk=kernel.dwk=0;}
                if(u_j<1) {kernel_main(u_j, hinv3_j, hinv4_j, &wk_j, &dwk_j, 1);} else {wk_j=dwk_j=0;}
                if(local.V_i<0 || isnan(local.V_i)) {local.V_i=0;}
                if(V_j<0 || isnan(V_j)) {V_j=0;}
                double sph_area = fabs(local.V_i*local.V_i*kernel.dwk + V_j*V_j*dwk_j); // effective face area //
                wk = 0.5 * (1 - 1/sqrt(1 + sph_area / (M_PI*kernel.r*kernel.r))); // corresponding geometric weight //

                if((wk <= 0)||(isnan(wk))) {continue;} // no point in going further, there's no physical weight here

                double wk_vec[AREA_WEIGHTED_SUM_ELEMENTS] = {0};
                wk_vec[0] = wk;
                if(kernel.dp[0]>0) {wk_vec[1]=wk*kernel.dp[0]/kernel.r; wk_vec[2]=0;} else {wk_vec[1]=0; wk_vec[2]=wk*kernel.dp[0]/kernel.r;}
                if(kernel.dp[1]>0) {wk_vec[3]=wk*kernel.dp[1]/kernel.r; wk_vec[4]=0;} else {wk_vec[3]=0; wk_vec[4]=wk*kernel.dp[1]/kernel.r;}
                if(kernel.dp[2]>0) {wk_vec[5]=wk*kernel.dp[2]/kernel.r; wk_vec[6]=0;} else {wk_vec[5]=0; wk_vec[6]=wk*kernel.dp[2]/kernel.r;}

                // ok worth initializing other variables we will use below
                #pragma omp atomic read
                InternalEnergy_j = SphP[j].InternalEnergy; // this can get modified below, so we need to read it thread-safe now
                for(k=0;k<3;k++) {
                    #pragma omp atomic read
                    Vel_j[k] = P[j].Vel[k]; // this can get modified below, so we need to read it thread-safe now
                }
                double InternalEnergy_j_0 = InternalEnergy_j, Mass_j_0 = Mass_j, rho_j_0 = rho_j, Vel_j_0[3]; for(k=0;k<3;k++) {Vel_j_0[k]=Vel_j[k];} // save initial values to use below
#ifdef METALS
                double Metallicity_j[NUM_METAL_SPECIES], Metallicity_j_0[NUM_METAL_SPECIES];
                for(k=0;k<NUM_METAL_SPECIES;k++) {
                    #pragma omp atomic read
                    Metallicity_j[k] = P[j].Metallicity[k]; // this can get modified below, so we need to read it thread-safe now
                    Metallicity_j_0[k] = Metallicity_j[k]; // save initial values to  use below
                }
#endif
                
                RsneKPC = RsneKPC_0;
                /* calculate cooling radius given density and metallicity in this annulus into which the ejecta propagate */
                if(loop_iteration < 2)
                {
#ifdef SN_KIMM_CEN_MODIFIED_MODEL
                    // density floor
                    double n0 = rho_j * density_to_n;
                    if (n0 < 0.001)
                    {
                        n0 = 0.001;
                    }

                    // metallicity floor
                    double z0 = Metallicity_j[0] / All.SolarAbundances[0];
                    if (z0 < 0.01)
                    {
                        z0 = 0.01;
                    }

                    // use scalings following Kimm & Cen (2014); Thornton et al. (1998)
                    const double n_dep = pow(n0, -4./17.);
                    const double z_dep = pow(z0, -0.28);
                    const double e_dep = pow(Esne51, 16./17.);
                    
                    // (Thornton+ 1998: p_terminal = 3.0e5 km/s per Msun [equiv. to Mcool = 895.5 Msun].)
                    // Injected radial momentum = \sqrt(2 * M_cool * E_SN), 
                    //Since the resoluton is increased, here I increase p_terminal = 5.0e5 km/s per Msun [equiv. to Mcool = 2487.5 Msun].
                    m_cooling = (895.5 / UNIT_MASS_IN_SOLAR) * (e_dep * n_dep * z_dep);
#else // NOT using SN_KIMM_CEN_MODIFIED_MODEL
                    double e0 = Esne51;
                    if(loop_iteration < 0) {e0=1;}
                    if(feedback_type_is_SNe == 1) {e0+=1;} else {e0=0.1;} // set to small number for non-SNe feedback
                    double n0 = rho_j*density_to_n;
                    if(n0 < 0.001) {n0=0.001;}
                    double z0 = Metallicity_j[0]/All.SolarAbundances[0];
                    if(z0 < 0.01) {z0 = 0.01;}
                    double z0_term=1.;
                    if(z0 < 1.) {z0_term = z0*sqrt(z0);} else {z0_term = z0;}
                    double nz_dep  = pow(n0 * z0_term , 0.14);
                    v_cooling = 210. * DMAX(nz_dep,0.5) / UNIT_VEL_IN_KMS;
                    m_cooling = 4.56e36 * e0 / (nz_dep*nz_dep * UNIT_MASS_IN_CGS);
#endif // SN_KIMM_CEN_MODIFIED_MODEL
                    RsneKPC = pow( 0.238732 * m_cooling/rho_j , 1./3. );
                }
                RsneKPC_3 = RsneKPC*RsneKPC*RsneKPC;
                // if loop_iteration==-1, this is a pre-calc loop to get the relevant weights for coupling //
                if(loop_iteration < 0)
                {
                    if(loop_iteration==-1) // the Area_weighted_sum quantities are computed on loop=-2; these quantities must be computed on loop=-1 (after Area_weighted_sums are computed)
                    {
                        /* calculate the corrected momentum vectors that we will actually use in the coupling proper */
                        double pnorm=0, pvec[3]={0}, vel_ba_2=0, cos_vel_ba_pcoupled=0;
                        for(k=0;k<3;k++)
                        {
                            double q = 0; int i1=2*k+1, i2=i1+1;
                            double q_i1 = fabs(local.Area_weighted_sum[i1]);
                            double q_i2 = fabs(local.Area_weighted_sum[i2]);
                            if((q_i1>MIN_REAL_NUMBER)&&(q_i2>MIN_REAL_NUMBER))
                            {
                                double rr = q_i2/q_i1;
                                double rr2 = rr * rr;
                                if(wk_vec[i1] != 0)
                                {
                                    // Eq. 10 of Hopkins+ 2018
                                    q += wk_norm * wk_vec[i1] * sqrt(0.5*(1.0+rr2));
                                } else {
                                    q += wk_norm * wk_vec[i2] * sqrt(0.5*(1.0+1.0/rr2));
                                }
                            } else {
                                q += wk_norm * (wk_vec[i1] + wk_vec[i2]);
                            }
                            pvec[k] = -q;
                            pnorm += pvec[k]*pvec[k];
                        }
                        pnorm = sqrt(pnorm);
                        /* now calculate the additional weights that are needed for energy terms */
                        for(k=0;k<3;k++)
                        {
                            double v_ba = (Vel_j[k] - local.Vel[k]) / All.cf_atime; // relative gas-star velocity //
                            vel_ba_2 += v_ba*v_ba; // magnitude of velocity vector (for corrected post-shock energies to distribute)
                            cos_vel_ba_pcoupled += v_ba * pvec[k]/pnorm; // direction of ejecta [after correction loop]
                        }
                        wk_vec[7] = pnorm * vel_ba_2; // phi_0 term : residual KE term from mass-coupling for {small, second-order} energy correction
                        wk_vec[8] = sqrt(pnorm * Mass_j) * cos_vel_ba_pcoupled; // beta_0 term : cross-term for momentum coupling effect on energy-coupling
                        wk_vec[9] = pnorm * cos_vel_ba_pcoupled / v_cooling; // calculate the beta term as if all particles hit terminal: more accurate result in that limit
                        wk_vec[10] = pnorm; // normalization (so that we can divide by its sum to properly normalize the beta_egy and beta_cool quantities)
                    }
                    for(k=0;k<AREA_WEIGHTED_SUM_ELEMENTS;k++) {out.Area_weighted_sum[k] += wk_vec[k];}
                    continue;
                }
                // NOW do the actual feedback calculation //
                wk *= wk_norm; // this way wk matches the value summed above for the weighting //

                if((wk <= 0)||(isnan(wk))) continue;

                /* define initial mass and ejecta velocity in this 'cone' */
                double pnorm = 0, pvec[3] = {0};
                for(k=0; k<3; k++)
                {
                    double q = 0; int i1=2*k+1, i2=i1+1;
                    double q_i1 = fabs(local.Area_weighted_sum[i1]);
                    double q_i2 = fabs(local.Area_weighted_sum[i2]);
                    if((q_i1>MIN_REAL_NUMBER)&&(q_i2>MIN_REAL_NUMBER))
                    {
                        double rr = q_i2/q_i1;
                        double rr2 = rr * rr;
                        if(wk_vec[i1] != 0)
                        {
                            q += wk_norm * wk_vec[i1] * sqrt(0.5*(1.0+rr2));
                        } else {
                            q += wk_norm * wk_vec[i2] * sqrt(0.5*(1.0+1.0/rr2));
                        }
                    } else {
                        q += wk_norm * (wk_vec[i1] + wk_vec[i2]);
                    }
                    pvec[k] = -q;
                    pnorm += pvec[k]*pvec[k];
                }
                pnorm = sqrt(pnorm); // this (vector norm) is the new 'weight function' for our purposes
                pnorm *= pnorm_sum; for(k=0;k<3;k++) {pvec[k] *= pnorm_sum;} // normalize following sum [10] to ensure faces sum to unity
                dM_ejecta_in = pnorm * local.Msne;
                double mj_preshock, massratio_ejecta;
                mj_preshock = Mass_j;
                massratio_ejecta = dM_ejecta_in / (dM_ejecta_in + Mass_j);

#ifdef SN_INJECTED_MOMENTUM_ACCOUNTING // Momentum/energy accounting
                // compute p_j[k] before injecting mass or momentum
                //   (in the frame of the star velocity v_i, i.e. where v_i = 0)
                double p_j_initial[3] = {0., 0., 0.};
                for (int k = 0; k < 3; ++k) {
                    p_j_initial[k] = Mass_j * (Vel_j[k] - local.Vel[k]);
                }
#endif // SN_INJECTED_MOMENTUM_ACCOUNTING

                /* inject actual mass from mass return */
                int couple_anything_but_scalar_mass_and_metals = 1; // key to indicate whether or not we actually need to do the next set of steps beyond pure scalar mass+metal couplings //
                if(P[j].Hsml<=0) {if(rho_j>0){rho_j*=(1+dM_ejecta_in/Mass_j);} else {rho_j=dM_ejecta_in*kernel.hinv3;}} else {rho_j+=kernel_zero*dM_ejecta_in*hinv3_j;}
                rho_j *= 1 + dM_ejecta_in/Mass_j; // inject mass at constant particle volume //
                Mass_j += dM_ejecta_in;
                out.M_coupled += dM_ejecta_in;
                
#ifdef METALS   /* inject metals */
                //for(k=0;k<NUM_METAL_SPECIES-NUM_AGE_TRACERS;k++) {Metallicity_j[k]=(1-massratio_ejecta)*Metallicity_j[k] + massratio_ejecta*local.yields[k];}
                for(k=0; k<NUM_METAL_SPECIES-NUM_AGE_TRACERS; k++) {
                    double old_metal = Metallicity_j[k];
                    Metallicity_j[k] = (1 - massratio_ejecta) * Metallicity_j[k] + massratio_ejecta * local.yields[k]/(local.Msne*1e10); // !!!!!!

                    // DEBUG: Track O16 injection
                    //if(k == 3 && local.yields[k] > 0) {
                    //    printf("SLUG_DEBUG: Metal injection to gas ID=%d: O16 from=%g to=%g, yield=%g, massratio=%g\n", 
                    //           P[j].ID, old_metal, Metallicity_j[k], local.yields[k]/(local.Msne*1e10), massratio_ejecta);
                    //    fflush(stdout);
                    //}
                }
#ifdef GALSF_FB_FIRE_AGE_TRACERS
                if(loop_iteration == 3) {for(k=NUM_METAL_SPECIES-NUM_AGE_TRACERS;k<NUM_METAL_SPECIES;k++) {Metallicity_j[k] += pnorm*local.yields[k]/Mass_j;}} // add age tracers in taking yields to mean MASS, so we can make it large without actually exchanging large masses
#ifndef GALSF_FB_FIRE_AGE_TRACERS_DISABLE_SURFACE_YIELDS
                if(loop_iteration != 3) {for(k=NUM_METAL_SPECIES-NUM_AGE_TRACERS;k<NUM_METAL_SPECIES;k++) {Metallicity_j[k]=(1-massratio_ejecta)*Metallicity_j[k] + massratio_ejecta*local.yields[k];}} // treat like any other yield when doing stellar mass exchange
#endif
#endif
#endif
                
                if(couple_anything_but_scalar_mass_and_metals)
                {
#ifdef SN_KIMM_CEN_MODIFIED_MODEL
                double wk_m_cooling = pnorm * m_cooling; // effective cooling mass for this particle

                // Kimm & Cen do not attempt to resolve the free-expansion phase, so no (+ 1) term inside sqrt(...)
                double boost_max = sqrt(wk_m_cooling / dM_ejecta_in); // terminal momentum boost-factor
                double boost_egycon = sqrt(mj_preshock / dM_ejecta_in); // energy-conserving limit for coupling through neighbors

                double mom_boost_fac = 1;
                if(feedback_type_is_SNe == 1)
                {
                    // appropriately re-weight boost to avoid energy conservation errors [cooling-limit]
                    boost_max *= psi_cool;
                    // appropriately re-weight boost to avoid energy conservation errors [energy-conserving-limit]
                    boost_egycon *= psi_egycon;

                    // limit to cooling case if egy-conserving exceeds terminal boost,
                    // or coupled mass short of cooling mass
                    if((wk_m_cooling < mj_preshock) || (boost_max < boost_egycon)) {
                        mom_boost_fac = boost_max;
                    } else {
                        mom_boost_fac = boost_egycon;
                    }
                }
#else // NOT using KIMM_CEN_MODIFIED_MODEL
#ifdef SLUG
#error "You must use KIMM_CEN_MODIFIED_MODEL with SLUG. Otherwise unphysical things will happen, since the ejecta energy/velocities are not tabulated."
#endif // SLUG
                /* inject momentum: account for ejecta being energy-conserving inside the cooling radius (or Hsml, if thats smaller) */
                double wk_m_cooling = pnorm * m_cooling; // effective cooling mass for this particle
                double boost_max = sqrt(1 + wk_m_cooling / dM_ejecta_in); // terminal momentum boost-factor
                double boost_egycon = sqrt(1 + mj_preshock / dM_ejecta_in); // energy-conserving limit for coupling through neighbors

                double mom_boost_fac = 1;
                if(feedback_type_is_SNe == 1)
                {
                    boost_max *= psi_cool; // appropriately re-weight boost to avoid energy conservation errors [cooling-limit]
                    boost_egycon *= psi_egycon; // appropriately re-weight boost to avoid energy conservation errors [energy-conserving-limit]
                    // limit to cooling case if egy-conserving exceeds terminal boost,
                    // or coupled mass short of cooling mass
                    double psi0 = 1; // factor to use below for velocity-limiter
                    if((wk_m_cooling < mj_preshock) || (boost_max < boost_egycon)) {
                        mom_boost_fac=boost_max;
                        psi0=DMAX(psi0,psi_cool);
                    } else {
                        mom_boost_fac=boost_egycon;
                        psi0=DMAX(psi0,psi_egycon);
                    }

                    if(mom_boost_fac < 1) {mom_boost_fac=1;} // impose lower limit of initial ejecta momentum

                    // finally account for simple physical limiter: if particle moving away faster than cooling terminal velocity, can't reach that velocity //
                    double vcool = DMIN(v_cooling/psi0 , v_ejecta_eff/mom_boost_fac); // effective velocity at stalling/cooling radius
                    double dv_dp_phys = 0; for(k=0;k<3;k++) {dv_dp_phys += (1-massratio_ejecta) * (kernel.dp[k]/kernel.r) * ((local.Vel[k] - Vel_j[k])/All.cf_atime);} // recession velocity of particle from SNe
                    double v_cooling_lim = DMAX( vcool , dv_dp_phys ); // cooling vel can't be smaller than actual vel (note: negative dvdp here automatically returns vcool, as desired)
                    double boostfac_max = DMIN(1000. , v_ejecta_eff/v_cooling_lim); // boost factor cant exceed velocity limiter - if recession vel large, limits boost
                    if(mom_boost_fac > boostfac_max) {mom_boost_fac = boostfac_max;} // apply limiter
                }
#endif //KIMM_CEN_MODIFIED_MODEL

                /* save summation values for outputs */
                dP = local.unit_mom_SNe / Mass_j * pnorm;
                dP_sum += dP; dP_boost_sum += dP * mom_boost_fac;

                /* actually do the injection */
		        // this gives the appropriately-normalized tap-able momentum from the energy-conserving solution
                double vel_prefactor =  mom_boost_fac * massratio_ejecta * (All.cf_atime*v_ejecta_eff) / pnorm;

#ifdef SN_MOMENTUM_LIMITER		
                // This appears to be necessary when SLUG is used for feedback, since it is
                //   possible (and happens frequently) to have multiple SN/particle/timestep.

                // Limit the maximum injected radial momentum to 6.0e5 Msun km/s
                //Since the injected momentum is increased by 5/3 times, we also increase the maximum here
                double prev_vel_prefactor = vel_prefactor * 1.0;
                const double max_radial_mom = 6.0e5 / (UNIT_VEL_IN_KMS * UNIT_MASS_IN_SOLAR);
                const double max_vel_prefactor = max_radial_mom / P[j].Mass;
                vel_prefactor = DMIN(vel_prefactor, max_vel_prefactor);
		
                // Limit delta_v to < 5% of speed of light
                // (in case something has gone badly wrong, such as injecting momentum into a particle with tiny mass)
                const double max_delta_v = 0.05 * C_LIGHT_CODE; // ~15,000 km/s
                vel_prefactor = DMIN(vel_prefactor, max_delta_v);
                
#endif // SN_MOMENTUM_LIMITER
		
                double KE_initial = 0, KE_final = 0;
                for(k=0; k<3; k++)
                {
                    double d_vel = vel_prefactor * pvec[k] + massratio_ejecta*(local.Vel[k] - Vel_j[k]); // local.Vel term from extra momentum of moving star, Vel_j term from going from momentum to velocity boost with added mass
                    KE_initial += Vel_j[k]*Vel_j[k]; Vel_j[k] += d_vel; KE_final += Vel_j[k]*Vel_j[k];
                }
		
                /* now calculate the residual energy and add it as thermal */
                KE_initial *= 0.5 * mj_preshock * All.cf_a2inv;
                KE_final *= 0.5 * Mass_j * All.cf_a2inv;
                double E_sne_initial = pnorm * Energy_injected_codeunits;
                double d_Egy_internal = KE_initial + E_sne_initial - KE_final;

#if !defined(SINGLE_STAR_FB_WINDS) /* (for single-star modules we ignore this b/c assume always trying to resolve R_cool) */
                //if(feedback_type_is_SNe == 1) /* if coupling radius > R_cooling, account for thermal energy loss in the post-shock medium: from Thornton et al. thermal energy scales as R^(-6.5) for R>R_cool. only use for SNe b/c scalings [like momentum] only apply there. over-cooling if code wants to do it will easily occur next timestep. */
                {
                    //if(d_Egy_internal < 0.5*E_sne_initial) {d_Egy_internal = 0.5*E_sne_initial;}
                    double r_eff_ij = kernel.r - Get_Particle_Size(j); /* get effective distance */
                    if(r_eff_ij > RsneKPC) {d_Egy_internal *= RsneKPC_3 / (r_eff_ij*r_eff_ij*r_eff_ij);} /* rescale the coupled energy as intended for the feedback mechanism */
                }
#endif
                d_Egy_internal /= Mass_j; // convert to specific internal energy, finally //
#ifndef MECHANICAL_FB_MOMENTUM_ONLY
                if(d_Egy_internal > 0)
                {
                    InternalEnergy_j += d_Egy_internal; E_coupled += d_Egy_internal;
                }
#endif
                    
                } // couple_anything_but_scalar_mass_and_metals
                    
                /* we updated variables that need to get assigned to element 'j' -- let's do it here in a thread-safe manner */
                #pragma omp atomic
                P[j].Mass += Mass_j - Mass_j_0; // finite mass update [delta difference added here, allowing for another element to update in the meantime]. done this way to ensure conservation.
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
                #pragma omp atomic
                SphP[j].MassTrue += Mass_j - Mass_j_0; // finite mass update
#endif
                if(rho_j_0 > 0) {
                    #pragma omp atomic
                    SphP[j].Density *= rho_j / rho_j_0; // inject mass at constant particle volume [no need to be exactly conservative here] //
                }
                for(k=0;k<3;k++) {
                    #pragma omp atomic
                    P[j].Vel[k] += Vel_j[k] - Vel_j_0[k]; // delta-update
                    #pragma omp atomic
                    SphP[j].VelPred[k] += Vel_j[k] - Vel_j_0[k]; // delta-update
                    #pragma omp atomic
                    P[j].dp[k] += Mass_j*Vel_j[k] - Mass_j_0*Vel_j_0[k]; // discrete momentum change
                }
                #pragma omp atomic
                SphP[j].InternalEnergy += InternalEnergy_j - InternalEnergy_j_0; // delta-update
                #pragma omp atomic
                SphP[j].InternalEnergyPred += InternalEnergy_j - InternalEnergy_j_0; // delta-update
                for(k=0;k<NUM_METAL_SPECIES;k++) {
                    #pragma omp atomic
                    P[j].Metallicity[k] += Metallicity_j[k] - Metallicity_j_0[k]; // delta-update
                }

#ifdef SN_INJECTED_MOMENTUM_ACCOUNTING // Momentum/energy accounting
                // compute p_j[k] after injecting mass or momentum
                //   (in the frame of the star velocity v_i, i.e. where v_i = 0)
                double p_j_final[3] = {0., 0., 0.};
                for (int k = 0; k < 3; ++k) {
                    p_j_final[k] = Mass_j * (Vel_j[k] - local.Vel[k]);
                }

                // compute dMomentum for particle j
                double dp_j[3] = {0., 0., 0.};
                for (int k = 0; k < 3; ++k) {
                    dp_j[k] = p_j_final[k] - p_j_initial[k];
                }

                // compute ||dMomentum|| == dmom_radial and add to cumulative total dMom
                double dp_j_normsq = 0.;
                for (int k = 0; k < 3; ++k) {
                    dp_j_normsq += dp_j[k]*dp_j[k];
                }
                const double dp_j_norm = sqrt(dp_j_normsq);

                // add to cumulative total dMom
                #pragma omp atomic
                out.injected_radial_momentum += dp_j_norm;
                #pragma omp atomic
                out.injected_thermal_energy += (d_Egy_internal * Mass_j);
                #pragma omp atomic
                out.pnorm += pnorm;
                for(int k=0; k<3; ++k) {
                    #pragma omp atomic
                    out.pvec[k] += pvec[k]; // sums to zero (convex hull closes, i.e. \sum_j A_j = 0)
                }
#endif // SN_INJECTED_MOMENTUM_ACCOUNTING
            } // for(n = 0; n < numngb; n++)
        } // while(startnode >= 0)
        if(mode == 1) {listindex++; if(listindex < NODELISTLENGTH) {startnode = DATAGET_NAME[target].NodeList[listindex]; if(startnode >= 0) {startnode = Nodes[startnode].u.d.nextnode;}}}    /* open it */
    } // while(startnode >= 0)

    /* Now collect the result at the right place */
    if(mode == 0) {out2particle_addFB(&out, target, 0, loop_iteration);} else {DATARESULT_NAME[target] = out;}
    return 0;
} // int addFB_evaluate

#endif // GALSF_USE_SNE_ONELOOP_SCHEME else


/* parent routine which calls the relevant loops */
void mechanical_fb_calc(int fb_loop_iteration)
{
#ifdef SN_INJECTED_MOMENTUM_ACCOUNTING
    for (int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
        P[i].SNe_InjectedMomentumThisStep = 0.;
        P[i].SNe_InjectedThermalEnergyThisStep = 0.;
        P[i].SNe_pnorm = 0.;
        for(int k = 0; k < 3; ++k) {
            P[i].SNe_pvec[k] = 0.;
        }
    }
#endif // SN_INJECTED_MOMENTUM_ACCOUNTING

    PRINT_STATUS(" ..mechanical feedback loop: iteration %d",fb_loop_iteration);
    #include "../system/code_block_xchange_perform_ops_malloc.h" /* this calls the large block of code which contains the memory allocations for the MPI/OPENMP/Pthreads parallelization block which must appear below */
    loop_iteration = fb_loop_iteration; /* sets the appropriate feedback type for the calls below */
    #include "../system/code_block_xchange_perform_ops.h" /* this calls the large block of code which actually contains all the loops, MPI/OPENMP/Pthreads parallelization */
    #include "../system/code_block_xchange_perform_ops_demalloc.h" /* this de-allocates the memory for the MPI/OPENMP/Pthreads parallelization block which must appear above */
    CPU_Step[CPU_SNIIHEATING] += measure_time(); /* collect timings and reset clock for next timing */

#ifdef DEBUG_RADIAL_MOMENTUM
    if (fb_loop_iteration == 0)
    {
        for (int i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
        {
            if (P[i].SNe_ThisTimeStep > 0)
            {
                std::ostringstream outputLine;
                outputLine << P[i].SNe_InjectedMomentumThisStep << " "
                           << P[i].SNe_InjectedThermalEnergyThisStep << " "
                           << P[i].SNe_ThisTimeStep << " "
                           << P[i].SNe_pnorm << " ";

                for(int k = 0; k < 3; ++k) {
                    outputLine << P[i].SNe_pvec[k] << " ";
                }

                spdlog::get("debug")->info(outputLine.str());
            }
        }
    }
#endif // DEBUG_RADIAL_MOMENTUM
}
#include "../system/code_block_xchange_finalize.h" /* de-define the relevant variables and macros to avoid compilation errors and memory leaks */


#endif /* GALSF_FB_MECHANICAL */
