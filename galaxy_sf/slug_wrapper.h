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

  void advanceToTime(double particle_age); // particle_age [yr]
  
  auto getNumberSNeThisTimestep() -> int;
  auto getNumberAliveStochasticStars() -> int;
  auto getYieldsThisTimestep() -> std::vector<double>;
  auto getBirthMass() -> double;
  auto getCurrentStellarMass() -> double;
  auto getPhotometryQH0() -> double; // ionising luminosity [photon/s]

  int numberSNeThisTimestep;
  std::vector<double> yieldsThisTimestep;

  // This is a pointer to the slug_cluster object
  slug_cluster *cluster;
};

#endif // SLUG_WRAPPER_H
