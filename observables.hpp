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
namespace EISH {
	/*
	 Store Spin state, Displacement state, Energy, Magnetization, ...
	 */
	struct SystemState {
		//constructor
		SystemState(int N);
    
    	void reset_sys_vars();
    	
    	void div_sys_vars(int N);
    
	    //structure elements
	    double sysE;
		double sysM;
		double sysMA;
		double sysMB;
		double sysSMA;
		double sysSMB;
		double sysDA;
		double sysDB;
		double sysSDA;
		double sysSDB;
		double sysMAMB;
		double sysDADB;
		double sysMADA;
		double sysMADB;
		double sysMBDA;
		double sysMBDB;
		double sysMAMBDA;
		double sysMAMBDB;
		double sysMADADB;
		double sysMBDADB;
		double sysMAMBDADB;
	    
	    
	    int pass;							// current pass
	
	    vector<int> SA;				// Spin state S in sublattice A
	    vector<int> SB;				// Spin state S in sublattice B
	    vector<double> DA;				// Distortion state in sublattice A
	    vector<double> DB;				// Distortion state in sublattice B
	};
}
}
#endif
