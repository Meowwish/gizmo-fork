#ifndef SLUG_WRAPPER_H
#define SLUG_WRAPPER_H

#include <slug.H>
#include <slug_cluster.H>
#include <slug_predefined.H>

constexpr auto do_stochastic_only = true;     // true is okay if only radioisotopes
constexpr auto minimum_stochastic_mass = 2.5; // minimum SNe mass is 9 Msun for Sukhbold16
constexpr auto stochastic_sampling_type = POISSON;
constexpr auto imf_type = "chabrier";
constexpr auto stellar_tracks = "modp020.dat";
constexpr auto spectral_synthesis = "sb99";
constexpr auto spectral_filter = "QH0";
constexpr auto yield_table = "SNII_Sukhbold16__AGB_Karakas16+Doherty14_nodecay";
constexpr auto recompute_yields_after_reconstruct = false; // *much* faster!

constexpr auto slug_cluster_internal_ID = 1;
constexpr auto slug_cluster_internal_time = 0.;

class slugWrapper
{

public:

  static slug_predefined *slug_globals;
  static int C12_index;
  static int C13_index;
  static int N14_index;
  static int N15_index;
  static int O16_index;
  static int O17_index;
  static int O18_index;
  static int Mg24_index;
  static int S32_index;
  static int Fe56_index;
  static int Ba138_index;
  static int Ce140_index;
  static int Eu153_index;

  // Constructor using particle_mass
  //  initialize a slug_cluster object with the given cluster mass (particle_mass)
  //  its internal ID is set to slug_cluster_internal_ID
  //  its internal time variable is set to slug_cluster_internal_time
  slugWrapper(double particle_mass)
      : yieldsThisTimestep_SNII(slug_globals->yields(yield_table)->get_niso()),
        yieldsThisTimestep_WR(slug_globals->yields(yield_table)->get_niso()),
        yieldsThisTimestep_AGB(slug_globals->yields(yield_table)->get_niso()),
        cluster(
            slug_cluster_internal_ID, particle_mass, slug_cluster_internal_time,
            slug_globals->imf(imf_type, minimum_stochastic_mass, stochastic_sampling_type),
            slug_globals->tracks(stellar_tracks),
            slug_globals->specsyn(spectral_synthesis, slug_globals->tracks(stellar_tracks),
                                slug_globals->imf(imf_type, minimum_stochastic_mass, stochastic_sampling_type)),
            slug_globals->filter_set(spectral_filter), nullptr, nullptr,
            slug_globals->yields(yield_table), nullptr,
            slug_globals->ostreams, nullptr, do_stochastic_only)
  {
  }

  // Method to reconstruct the slug_cluster object from a serialized buffer
  slugWrapper(slug_cluster_state_noyields &state)
      : yieldsThisTimestep_SNII(slug_globals->yields(yield_table)->get_niso()),
        yieldsThisTimestep_WR(slug_globals->yields(yield_table)->get_niso()),
        yieldsThisTimestep_AGB(slug_globals->yields(yield_table)->get_niso()),
        cluster(state,
                     slug_globals->imf(imf_type, minimum_stochastic_mass, stochastic_sampling_type),
                     slug_globals->tracks(stellar_tracks),
                     slug_globals->specsyn(spectral_synthesis, slug_globals->tracks(stellar_tracks),
                                         slug_globals->imf(imf_type, minimum_stochastic_mass, stochastic_sampling_type)),
                     slug_globals->filter_set(spectral_filter), nullptr, nullptr,
                     slug_globals->yields(yield_table), nullptr,
                     slug_globals->ostreams, nullptr, do_stochastic_only, recompute_yields_after_reconstruct)
  {
  }

  // method to save the slug_cluster object
  void serializeCluster(slug_cluster_state_noyields &state);

  // method to advance cluster object in time
  void advanceToTime(double particle_age); // particle_age [yr]
  
  //method to return tracks->max_age()
  double getMaxAge() const{
  return cluster.getMaxAge();
  }


  // accessor functions
  auto getNumberSNeThisTimestep() -> int;
  auto getNumberAliveStochasticStars() -> int;
  auto getYieldsThisTimestep_SNII() -> std::vector<double>;
  auto getYieldsThisTimestep_WR() -> std::vector<double>;
  auto getYieldsThisTimestep_AGB() -> std::vector<double>;
  auto getEjectaMassThisTimestep() -> double;
  auto getPhotometryQH0() -> double; // ionising luminosity [photon/s]

  // member variables
  int numberSNeThisTimestep;
  std::vector<double> yieldsThisTimestep_SNII;
  std::vector<double> yieldsThisTimestep_WR;
  std::vector<double> yieldsThisTimestep_AGB;
  slug_cluster cluster;
};

#endif // SLUG_WRAPPER_H
