#include "slug_feedback.hpp"

void slugComputeSNFeedback(int i)
{
    // use SLUG to determine whether a SN event has occured in the last timestep
    if (P[i].slug_state_initialized)
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

#ifdef SLUG_YIELDS
        // WARNING: implementation not complete!
        // (TODO: need to extract yields for specified isotopes into GIZMO arrays)
        auto yields = mySlugObject.getYieldsThisTimestep(); // solar mass
        assert(yields.size() == NUM_METAL_SPECIES);

        for (size_t j = 0; j < yields.size(); ++j)
        {
            P[i].Yields_ThisTimestep[j] = yields[j];
        }
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