#ifndef SLUG_WRAPPER_H
#define SLUG_WRAPPER_H

#include <slug.H>
#include <slug_cluster.H>
#include <slug_predefined.H>

class slugWrapper {
public:
  // Constructor
  slugWrapper() { cluster = nullptr; }

  // Destructor; this frees the slug_cluster object
  virtual ~slugWrapper() {
    delete cluster;
    cluster = nullptr;
  }

  // Method to construct the slug_cluster object from particle mass
  void constructCluster(double particle_mass);
  void reconstructCluster(slug_cluster_state_noyields &state);
  void serializeCluster(slug_cluster_state_noyields &state);

  auto advanceToTime(double particle_age) -> std::vector<double>; // particle_age [yr]
  
  auto getStochasticSN() -> int;
  auto getBirthMass() -> double;
  auto getStellarMass() -> double;
  auto getPhotometryQH0() -> double; // ionising luminosity [photon/s]

  // This is a pointer to the slug_cluster object
  slug_cluster *cluster;
};

#endif // SLUG_WRAPPER_H
