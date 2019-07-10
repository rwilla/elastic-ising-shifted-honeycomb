/**
 * @file observables.hpp
 *
 * Classes related to definition of the state observables
 */

#ifndef OBSERVABLES_HPP
#define OBSERVABLES_HPP

#include <iostream>
#include <vector>

using namespace std;

namespace observables {

/*
 Store Spin state, Displacement state, Energy, Magnetization, ...
 */
struct SystemState {
	//constructor
	SystemState(int N);

	//print functions
	void print_E();

    
    //structure elements
    double E;							// Energy per spin

    double MA;							// Magnetization of sublattice A
    double MB;							// Magnetization of sublattice B

    double sMA;							// Staggered Magnetization A
    double sMB;							// Staggered Magnetization B

    double Lnem;						// nematic order MA * MB
    double Ldel;						// DA - DB
    
    int pass;							// current pass

    vector<int> SA;				// Spin state S in sublattice A
    vector<int> SB;				// Spin state S in sublattice B
    vector<float> DA;				// Distortion state in sublattice A
    vector<float> DB;				// Distortion state in sublattice B
};

}

#endif
