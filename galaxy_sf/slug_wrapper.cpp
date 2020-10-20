#include "slug_wrapper.h"

slug_predefined *slugWrapper::slug_globals = nullptr;

void slugWrapper::serializeCluster(slug_cluster_state_noyields &state)
{
  cluster.serializeToStructWithoutYields(state);
}

void slugWrapper::advanceToTime(double particle_age)
{
  const std::vector<double> yields_t0 = cluster.get_yield();
  const int numberSNe_t0 = cluster.get_stoch_sn();

  cluster.advance(particle_age);

  const std::vector<double> yields_t1 = cluster.get_yield();
  const int numberSNe_t1 = cluster.get_stoch_sn();

  for (size_t i = 0; i < yieldsThisTimestep.size(); ++i)
  {
    yieldsThisTimestep[i] = std::max(yields_t1[i] - yields_t0[i], 0.0);
  }

  numberSNeThisTimestep = numberSNe_t1 - numberSNe_t0;
}

auto slugWrapper::getNumberSNeThisTimestep() -> int
{
  // get the number of supernovae that went off during this timestep
  // (this is computed within advanceToTime(double time))
  return numberSNeThisTimestep;
}

auto slugWrapper::getYieldsThisTimestep() -> std::vector<double>
{
  return yieldsThisTimestep;
}

auto slugWrapper::getEjectaMassThisTimestep() -> double
{
  double ejectaMass = 0.;
  // we sum the yields of all elements to compute the total ejecta mass
  for (size_t i = 0; i < yieldsThisTimestep.size(); ++i)
  {
    ejectaMass += yieldsThisTimestep[i];
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