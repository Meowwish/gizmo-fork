#include "slug_wrapper.h"

slug_predefined *slugWrapper::slug_globals = nullptr;

void slugWrapper::serializeCluster(slug_cluster_state_noyields &state)
{
  cluster.serializeToStructWithoutYields(state);
}

void slugWrapper::advanceToTime(double particle_age)
{
  const std::vector<double> yields_SNII_t0 = cluster.get_yield(YIELDS_SNII);
  const std::vector<double> yields_WR_t0 = cluster.get_yield(YIELDS_MASSIVE_WIND);
  const std::vector<double> yields_AGB_t0 = cluster.get_yield(YIELDS_AGB);
  const int numberSNe_t0 = cluster.get_stoch_sn();
  
  
  //printf("SLUG AGB DEBUG: current AGB yield size: %zu\n", yields_AGB_t0.size());
  //bool all_zeros_AGB_t0 = true;
  //double sum_AGB_t0 = 0.0;
  //for (size_t i = 0; i < yields_AGB_t0.size(); ++i) {
  //  if (yields_AGB_t0[i] != 0.0) all_zeros_AGB_t0 = false;
  //  sum_AGB_t0 += yields_AGB_t0[i];
  //}
  //printf("SLUG AGB DEBUG: current AGB yields all zeros? %s, Sum: %g\n", 
  //        all_zeros_AGB_t0 ? "Yes" : "No", sum_AGB_t0);
  
  //printf("SLUG WR DEBUG: current WR yield size: %zu\n", yields_WR_t0.size());
  //bool all_zeros_WR_t0 = true;
  //double sum_WR_t0 = 0.0;
  //for (size_t i = 0; i < yields_WR_t0.size(); ++i) {
  //  if (yields_WR_t0[i] != 0.0) all_zeros_WR_t0 = false;
  //  sum_WR_t0 += yields_WR_t0[i];
  //}
  //printf("SLUG WR DEBUG: current WR yields all zeros? %s, Sum: %g\n", 
  //        all_zeros_WR_t0 ? "Yes" : "No", sum_WR_t0);
  
  //printf("SLUG DEBUG (slug.wrapper.cpp): beforeAdvance yields_t0[9](C12) = %g, yields_t0[14](O16) = %g\n", 
  //       yields_t0[9], yields_t0[14]);
  
  cluster.advance(particle_age);
  
  const std::vector<double> yields_SNII_t1 = cluster.get_yield(YIELDS_SNII);
  const std::vector<double> yields_WR_t1 = cluster.get_yield(YIELDS_MASSIVE_WIND);
  const std::vector<double> yields_AGB_t1 = cluster.get_yield(YIELDS_AGB);

  const int numberSNe_t1 = cluster.get_stoch_sn();
  
  //if (numberSNe_t1 - numberSNe_t0 > 0) {
  //  printf("SLUG DEBUG (slug.wrapper.cpp): afterAdvance yields_t1[9](C12) = %.2f, yields_t1[14](O16) = %.2f\n", 
  //       yields_t1[9]*1e10, yields_t1[14]*1e10);
  //  printf("SLUG DEBUG (slug.wrapper.cpp): SNe_A = %.1f, SNe_B = %.1f\n", numberSNe_t0, numberSNe_t1);}
  
  for (size_t i = 0; i < yieldsThisTimestep_SNII.size(); ++i)
  {
    yieldsThisTimestep_SNII[i] = std::max(yields_SNII_t1[i] - yields_SNII_t0[i], 0.0);
    yieldsThisTimestep_WR[i] = std::max(yields_WR_t1[i] - yields_WR_t0[i], 0.0);
    yieldsThisTimestep_AGB[i] = std::max(yields_AGB_t1[i] - yields_AGB_t0[i], 0.0);
    //printf("SLUG AGB DEBUG: yields before advanced: %g, after advanced: %g, difference = %g\n",yields_AGB_t0[i],yields_AGB_t1[i],yields_AGB_t1[i] - yields_AGB_t0[i]);
  }
  
  //printf("SLUG AGB DEBUG: advanced AGB yield size: %zu\n", yields_AGB_t1.size());
  //bool all_zeros_AGB_t1 = true;
  //double sum_AGB_t1 = 0.0;
  //for (size_t i = 0; i < yields_AGB_t1.size(); ++i) {
  //  if (yields_AGB_t1[i] != 0.0) all_zeros_AGB_t1 = false;
  //  sum_AGB_t1 += yields_AGB_t1[i];
  //}
  //printf("SLUG AGB DEBUG: advanced AGB yields all zeros? %s, Sum: %g\n", 
  //        all_zeros_AGB_t1 ? "Yes" : "No", sum_AGB_t1);

  //printf("SLUG WR DEBUG: current WR yield size: %zu\n", yields_WR_t1.size());
  //bool all_zeros_WR_t1 = true;
  //double sum_WR_t1 = 0.0;
  //for (size_t i = 0; i < yields_WR_t0.size(); ++i) {
  //  if (yields_WR_t1[i] != 0.0) all_zeros_WR_t1 = false;
  //  sum_WR_t1 += yields_WR_t1[i];
  //}
  //printf("SLUG WR DEBUG: current WR yields all zeros? %s, Sum: %g\n", 
  //        all_zeros_WR_t1 ? "Yes" : "No", sum_WR_t1);
  
  //if (numberSNe_t1 - numberSNe_t0 > 0) {
  //  printf("SLUG DEBUG (slug.wrapper.cpp): yieldsThisTimestep[9](C12) = %.2f, yieldsThisTimestep[14](O16) = %.2f\n",
  //       yieldsThisTimestep[9]*1e10, yieldsThisTimestep[14]*1e10);}
  
  numberSNeThisTimestep = numberSNe_t1 - numberSNe_t0;
}

auto slugWrapper::getNumberSNeThisTimestep() -> int
{
  // get the number of supernovae that went off during this timestep
  // (this is computed within advanceToTime(double time))
  return numberSNeThisTimestep;
}

auto slugWrapper::getYieldsThisTimestep_SNII() -> std::vector<double>
{
  return yieldsThisTimestep_SNII;
}
auto slugWrapper::getYieldsThisTimestep_WR() -> std::vector<double>
{
  return yieldsThisTimestep_WR;
}
auto slugWrapper::getYieldsThisTimestep_AGB() -> std::vector<double>
{
  return yieldsThisTimestep_AGB;
}

auto slugWrapper::getEjectaMassThisTimestep() -> double
{
  double ejectaMass = 0.;
  // we sum the yields of all elements to compute the total ejecta mass
  for (size_t i = 0; i < yieldsThisTimestep_SNII.size(); ++i)
  {
    ejectaMass += yieldsThisTimestep_SNII[i];
    ejectaMass += yieldsThisTimestep_WR[i];
    ejectaMass += yieldsThisTimestep_AGB[i];
  }

  return ejectaMass; // solar masses
}

auto slugWrapper::getNumberAliveStochasticStars() -> int
{
  // get the number of stars (above the minimum_stochastic_mass mass threshold) that are still alive
  return cluster.get_nstars();
}

auto slugWrapper::getPhotometryQH0() -> double
{
  // returns the number of photons with energy > 13.6 eV
  //   produced at the current time by stars in the cluster
  return cluster.get_photometry()[0];
}