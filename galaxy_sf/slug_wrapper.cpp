#include "slug_wrapper.h"

constexpr auto do_stochastic_only = true; // true is okay if only radioisotopes
constexpr auto minimum_stochastic_mass = 9.0; // minimum SNe mass is 9 Msun for Sukhbold16
constexpr auto stochastic_sampling_type = POISSON;
constexpr auto imf_type = "chabrier";
constexpr auto stellar_tracks = "modp020.dat";
constexpr auto spectral_synthesis = "sb99";
constexpr auto spectral_filter = "QH0";
constexpr auto yield_table = "SNII_Sukhbold16_nodecay";
constexpr auto compute_yields = true;

constexpr auto slug_cluster_internal_ID = 1;
constexpr auto slug_cluster_internal_time = 0.;

void slugWrapper::constructCluster(double particle_mass)
{
  // initialize a slug_cluster object with the given cluster mass (particle_mass)
  // its internal ID is set to slug_cluster_internal_ID
  // its internal time variable is set to slug_cluster_internal_time
  cluster = new slug_cluster(
      slug_cluster_internal_ID, particle_mass, slug_cluster_internal_time,
      slug_predef.imf(imf_type, minimum_stochastic_mass, stochastic_sampling_type),
      slug_predef.tracks(stellar_tracks),
      slug_predef.specsyn(spectral_synthesis, slug_predef.tracks(stellar_tracks),
                          slug_predef.imf(imf_type, minimum_stochastic_mass, stochastic_sampling_type)),
      slug_predef.filter_set(spectral_filter), nullptr, nullptr,
      slug_predef.yields(yield_table), nullptr,
      slug_predef.ostreams, nullptr, do_stochastic_only, compute_yields);
}

// Method to reconstruct the slug_cluster object from a serialized buffer
void slugWrapper::reconstructCluster(slug_cluster_state_noyields &state)
{
  // reconstruct a slug_cluster from a slug_cluster_state_noyields struct
  cluster = new slug_cluster(
      state,
      slug_predef.imf(imf_type, minimum_stochastic_mass, stochastic_sampling_type),
      slug_predef.tracks(stellar_tracks),
      slug_predef.specsyn(spectral_synthesis, slug_predef.tracks(stellar_tracks),
                          slug_predef.imf(imf_type, minimum_stochastic_mass, stochastic_sampling_type)),
      slug_predef.filter_set(spectral_filter), nullptr, nullptr,
      slug_predef.yields(yield_table), nullptr,
      slug_predef.ostreams, nullptr, do_stochastic_only, compute_yields);
}

void slugWrapper::serializeCluster(slug_cluster_state_noyields &state)
{
  cluster->serializeToStructWithoutYields(state);
}

void slugWrapper::advanceToTime(double particle_age)
{
  const std::vector<double> yields_t0 = cluster->get_yield();
  const int numberSNe_t0 = cluster->get_stoch_sn();

  cluster->advance(particle_age);

  const std::vector<double> yields_t1 = cluster->get_yield();
  const int numberSNe_t1 = cluster->get_stoch_sn();

  // TODO: yieldsThisTimestep needs to be initialized to the size of the yield table in use
  for(size_t i=0; i < yieldsThisTimestep.size(); ++i) {
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

auto slugWrapper::getNumberAliveStochasticStars() -> int
{
  // get the number of stars (above the minimum_stochastic_mass mass threshold) that are still alive
  return cluster->get_nstars();
}

auto slugWrapper::getBirthMass() -> double
{
  // get the zero-age main sequence stellar mass for all stars in the cluster
  return cluster->get_birth_mass();
}

auto slugWrapper::getCurrentStellarMass() -> double
{
  // this is not identical to the birth mass due to mass loss in stellar winds
  return cluster->get_stellar_mass();
}

auto slugWrapper::getPhotometryQH0() -> double
{
  // returns the number of photons with energy > 13.6 eV
  //   produced at the current time by stars in the cluster
  return cluster->get_photometry()[0];
}