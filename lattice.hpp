/**
 * @file lattice.hpp
 *
 * Classes related to definition of the lattice
 */

#ifndef LATTICE_HPP
#define LATTICE_HPP

#include <vector>

using namespace std;

namespace lattice {

/*
 Store indices to identify nearest neighbors.
 */
struct NearestNeighbors  {
    
    NearestNeighbors(int N);
    
    vector<int> AB1;					// nearest neighbor 1 to a spin A belonging to sublattice B in layer below
    vector<int> AB2;					// nearest neighbor 2 to a spin A belonging to sublattice B in layer below
    vector<int> AB3;					// nearest neighbor 3 to a spin A belonging to sublattice B in layer below
    vector<int> AB4;					// nearest neighbor 1 to a spin A belonging to sublattice B in layer above
    vector<int> AB5;					// nearest neighbor 2 to a spin A belonging to sublattice B in layer above
    vector<int> AB6;					// nearest neighbor 3 to a spin A belonging to sublattice B in layer above

    vector<int> AA1;					// nearest neighbor to a spin A spin belonging to sublattice A below
    vector<int> AA2;					// nearest neighbor to a spin A spin belonging to sublattice A above

    vector<int> BA1;					// nearest neighbor 1 to a spin B belonging to sublattice A in layer below
    vector<int> BA2;					// nearest neighbor 2 to a spin B belonging to sublattice A in layer below
    vector<int> BA3;					// nearest neighbor 3 to a spin B belonging to sublattice A in layer below
    vector<int> BA4;					// nearest neighbor 1 to a spin B belonging to sublattice A in layer above
    vector<int> BA5;					// nearest neighbor 2 to a spin B belonging to sublattice A in layer above
    vector<int> BA6;					// nearest neighbor 3 to a spin B belonging to sublattice A in layer above

    vector<int> BB1;					// nearest neighbor to a spin B spin belonging to sublattice B below
    vector<int> BB2;					// nearest neighbor to a spin B spin belonging to sublattice B above

	vector<int>sweeporder;				// order in which spins are visited
};

}

#endif
