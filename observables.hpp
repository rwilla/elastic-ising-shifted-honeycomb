/**
 * @file observables.hpp
 *
 * Classes related to definition of the state observables
 */

#ifndef OBSERVABLES_HPP
#define OBSERVABLES_HPP

#include <iostream>
#include <vector>

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

    std::vector<int> SA;				// Spin state S in sublattice A
    std::vector<int> SB;				// Spin state S in sublattice B
    std::vector<float> DA;				// Distortion state in sublattice A
    std::vector<float> DB;				// Distortion state in sublattice B
};

}

#endif
