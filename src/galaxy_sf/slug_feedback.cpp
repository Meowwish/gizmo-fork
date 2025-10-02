#include "slug_feedback.hpp"

void slugComputeSNFeedback(int i)
{
    // use SLUG to determine whether a SN event has occured in the last timestep
    if (P[i].StellarAge > 0)
    {
        // create slug object
        slugWrapper mySlugObject(P[i].slug_state);

        // advance slug object in time
        // [the slug object should NOT be advanced in time anywhere else in the code,
        //     otherwise the yields and SNe events will not be accounted for.]
        double cluster_age_in_years = (All.Time - P[i].StellarAge) * UNIT_TIME_IN_YR;
        mySlugObject.advanceToTime(cluster_age_in_years);

        P[i].SNe_ThisTimeStep = mySlugObject.getNumberSNeThisTimestep();         // dimensionless
        P[i].EjectaMass_ThisTimestep = mySlugObject.getEjectaMassThisTimestep(); // solar mass

        // keep track of the cumulative number of SNe produced by this particle
        P[i].SNe_Cumulative += P[i].SNe_ThisTimeStep;
        //if (P[i].SNe_ThisTimeStep>0){printf("SNe = %.1f, Cumulative = %.1f\n", P[i].SNe_ThisTimeStep, P[i].SNe_Cumulative);}

#ifdef SLUG_YIELDS
        // WARNING: implementation not complete!
        // (TODO: need to extract yields for specified isotopes into GIZMO arrays)
        auto yields_SNII = mySlugObject.getYieldsThisTimestep_SNII(); // solar mass
        auto yields_WR = mySlugObject.getYieldsThisTimestep_WR(); // solar mass
        auto yields_AGB = mySlugObject.getYieldsThisTimestep_AGB(); // solar mass
        //if (P[i].SNe_ThisTimeStep>0){printf("SLUG DEBUG: yields.size() = %zu\n", yields_AGB.size());}

        int C12_index = slugWrapper::slug_globals->yields(yield_table)->isotope_index(6,12);
        int C13_index = slugWrapper::slug_globals->yields(yield_table)->isotope_index(6,13);
        int N14_index = slugWrapper::slug_globals->yields(yield_table)->isotope_index(7,14);
        int N15_index = slugWrapper::slug_globals->yields(yield_table)->isotope_index(7,15);
        int O16_index = slugWrapper::slug_globals->yields(yield_table)->isotope_index(8,16);
        int O17_index = slugWrapper::slug_globals->yields(yield_table)->isotope_index(8,17);
        int O18_index = slugWrapper::slug_globals->yields(yield_table)->isotope_index(8,18);
        int Mg24_index = slugWrapper::slug_globals->yields(yield_table)->isotope_index(12,24);
        int S32_index = slugWrapper::slug_globals->yields(yield_table)->isotope_index(16,32);
        int Fe56_index = slugWrapper::slug_globals->yields(yield_table)->isotope_index(26,56);
        int Ba138_index = slugWrapper::slug_globals->yields(yield_table)->isotope_index(56,138);
        int Ce140_index = slugWrapper::slug_globals->yields(yield_table)->isotope_index(58,140);
        int Eu153_index = slugWrapper::slug_globals->yields(yield_table)->isotope_index(63,153);
        //printf("SLUG DEBUG: C12_index = %d, N14_index = %d, O16_index = %d\n", 
        //    C12_index, N14_index, O16_index);

        std::vector<int> indices_SNII = {C12_index, C13_index, N14_index, N15_index, O16_index, O17_index, O18_index, Mg24_index, S32_index, Fe56_index, Ba138_index, Ce140_index, Eu153_index};
        std::vector<int> indices_WR = {C12_index, C13_index, N14_index, N15_index, O16_index, O17_index, O18_index, Mg24_index, S32_index, Fe56_index, Ba138_index, Ce140_index, Eu153_index};
        std::vector<int> indices_AGB = {C12_index, C13_index, N14_index, N15_index, O16_index, O17_index, O18_index, Mg24_index, S32_index, Fe56_index, Ba138_index, Ce140_index, Eu153_index};
        for (size_t idx = 0; idx < indices_SNII.size(); ++idx)
        {
            size_t isotopeIndex = indices_SNII[idx];
            //size_t j = idx + 1; // Adjusting idx to start from 1 instead of 0
            size_t j = idx + 1; // Adjusting idx to start from 1 instead of 0
            //printf("isotopeIndex[%zu] = %zu, \\\n", idx, isotopeIndex);

            if (isotopeIndex < yields_SNII.size()) {
                P[i].Yields_ThisTimestep[j] = yields_SNII[isotopeIndex]; // Use isotopeIndex to access the correct yield
                //if (j == 3 and P[i].SNe_ThisTimeStep>0){printf("SNe = %.1f, Yields_ThisTimestep_SNII[%zu](O16) = %g\n", P[i].SNe_ThisTimeStep, j, P[i].Yields_ThisTimestep[j]);}
            } else {
                std::cerr << "Index " << isotopeIndex << " is out of bounds for the SNII yields vector." << std::endl;
            }
        }
        for (size_t idx = 0; idx < indices_WR.size(); ++idx)
        {
            size_t isotopeIndex = indices_WR[idx];
            size_t j = idx + 14; // Adjusting idx to start from 6 instead of 0
            //printf("isotopeIndex[%zu] = %zu, \\\n", idx, isotopeIndex);

            if (isotopeIndex < yields_WR.size()) {
                P[i].Yields_ThisTimestep[j] = yields_WR[isotopeIndex]; // Use isotopeIndex to access the correct yield
                //if (j == 3 + 6 and P[i].Yields_ThisTimestep[j] > 0){printf("SNe = %.1f, Yields_ThisTimestep_WR[%zu](O16) = %g\n", P[i].SNe_ThisTimeStep, j, P[i].Yields_ThisTimestep[j]);}
            } else {
                std::cerr << "Index " << isotopeIndex << " is out of bounds for the WR yields vector." << std::endl;
            }
        }
        for (size_t idx = 0; idx < indices_AGB.size(); ++idx)
        {
            size_t isotopeIndex = indices_AGB[idx];
            size_t j = idx + 27; // Adjusting idx to start from 6 instead of 0
            //printf("isotopeIndex[%zu] = %zu, \\\n", idx, isotopeIndex);

            if (isotopeIndex < yields_AGB.size()) {
                P[i].Yields_ThisTimestep[j] = yields_AGB[isotopeIndex]; // Use isotopeIndex to access the correct yield
                // When AGB iron yield (assumed at indices_AGB index 5, j==24) is > 0, output debug info.
                /*
                if (idx == 5 && P[i].Yields_ThisTimestep[j] > 0) {
                    double cluster_age_in_years = (All.Time - P[i].StellarAge) * UNIT_TIME_IN_YR;
                    std::cout << "DEBUG: Particle " << i << " information:" << std::endl;
                    std::cout << "  StellarAge: " << P[i].StellarAge << std::endl;
                    std::cout << "  SNe_ThisTimeStep: " << P[i].SNe_ThisTimeStep << std::endl;
                    std::cout << "  EjectaMass_ThisTimestep: " << P[i].EjectaMass_ThisTimestep << std::endl;
                    std::cout << "  Yields_ThisTimestep: ";
                    std::cout << P[i].Yields_ThisTimestep[j] << std::endl;
                    std::cout << "  Current Time: " << All.Time 
                            << ", Advanced Time: " << cluster_age_in_years << " years" << std::endl;
                    std::exit(0); // Immediately exit the simulation program
                }
                */
            } else {
                std::cerr << "Index " << isotopeIndex << " is out of bounds for the AGB yields vector." << std::endl;
            }
        }
        //assert(yields.size() == NUM_METAL_SPECIES);

        /*for (size_t j = 0; j < yields.size(); ++j)
        {
            P[i].Yields_ThisTimestep[j] = yields[j];
        }*/ // Old
#endif // SLUG_YIELDS

        // serialize slug object
        mySlugObject.serializeCluster(P[i].slug_state);

        // check whether all stochastic stars have died
        if (mySlugObject.getNumberAliveStochasticStars() == 0)
        {
            // if so, mark the object as inactive
            P[i].slug_state_initialized = false;
        }
    } // mySlugObject deallocated automatically
}

auto slugComputeIonizingPhotons(int i) -> double
{
    // compute number of ionizing photons via SLUG
    double N_photons = 0.; // units == [s^-1]

    if (P[i].slug_state_initialized)
    {
        // re-create slug object
        slugWrapper mySlugObject(P[i].slug_state);

        // advance to current time (*but do NOT save, otherwise we will lose SN events!*)
        double cluster_age_in_years = (All.Time - P[i].StellarAge) * UNIT_TIME_IN_YR;
        mySlugObject.advanceToTime(cluster_age_in_years);

        // compute ionizing photons
        N_photons = mySlugObject.getPhotometryQH0();
    }
    return N_photons;
}