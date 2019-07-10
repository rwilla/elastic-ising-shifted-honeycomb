/**
 * @file iparameters.hpp
 *
 * Defines default initial parameters
 * and reads the input line for changing these parameters
 */

#ifndef INITIALPARAMETERS_HPP
#define INITIALPARAMETERS_HPP


#include <string>

using namespace std;

namespace iparameters {

/*
 Store indices to identify nearest neighbors.
 */
struct InitialParameters  {
    
    InitialParameters();
    
    
    //system
    int Lx;				// in-plane size of the cube lattice
	int Ly;				// in-plane size of the cube lattice
	int Lz;				// oo-plane size of the cube lattice
	int N;				// Number of spins
	int rs;						// random seed

	//physical parameters
	double T; 			// Temperature
	double Hx;			// Field strength along z
	double Hy;			// Field strength along z
	double Hz;			// Field strength along z
	
	//spin interactions
	double Jcc;			// in-plane ferromagnetic coupling (antiferromagnetic if Jxy < 0)
	double Jab;			// oo-plane ferromagnetic coupling (antiferromagnetic if Joo < 0)
	
	double K;			// c-axis anisotropy
	double L;			// in-plane anisotropy
	
	//elastic interactions
	double eps;						
	
	//simulation paramters
	int passes;					// Default # of passes
	string printC_Q;			// flag defining if the correlation function shall be printed
	string printLastConf_Q;		// flag defining if the last configuration shall be printed
	
};

}

#endif
