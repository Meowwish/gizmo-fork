#include "slug_sfr.hpp"

void slugFormStar(int i)
{
    const double cluster_mass = P[i].Mass * UNIT_MASS_IN_SOLAR;
    slugWrapper mySlugObject(cluster_mass);
    mySlugObject.serializeCluster(P[i].slug_state);
    P[i].slug_state_initialized = true;
}