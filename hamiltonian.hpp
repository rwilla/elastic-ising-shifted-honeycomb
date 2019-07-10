/**
 * @file hamiltonian.hpp
 *
 * Defines the contributions to the total hamiltonian
 * different lattice realizations may be implemented here. Therefore each hamiltonian
 * is prepend by an identified, e.g. EISH elastic ising shifted honeycomb
 * each spin hamiltonian gets its own namespace
 */

#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP


#include <string>
#include "lattice.hpp"
#include "observables.hpp"
#include "iparameters.hpp"

using namespace std;

namespace hamiltonian {

// define elastic ising hamiltonian for the shifted honeycomb lattice
// Ising spin with 6 neighbors on other sublattice and 2 neighbors on the same sublattice
// EISH = elastic ising shifted honeycomb
namespace EISH {
//onsite contribution
double onsite(const iparameters::InitialParameters& ip, double S, double D, int spin_index);

//Full hamiltonian for a spin on sublattice A
double siteA(const iparameters::InitialParameters& ip, double S, double D, int spin_index, const lattice::NearestNeighbors& nn, observables::SystemState& s);

//Full hamiltonian for a spin on sublattice B
double siteB(const iparameters::InitialParameters& ip, double S, double D, int spin_index, const lattice::NearestNeighbors& nn, observables::SystemState& s);

//Nearest neighbor interaction for a spin on sublattice A
double nnA(const iparameters::InitialParameters& ip, double S, double D, int spin_index, const lattice::NearestNeighbors& nn, observables::SystemState& s);

//Nearest neighbor interaction for a spin on sublattice B
double nnB(const iparameters::InitialParameters& ip, double S, double D, int spin_index, const lattice::NearestNeighbors& nn, observables::SystemState& s);

}

}

#endif
