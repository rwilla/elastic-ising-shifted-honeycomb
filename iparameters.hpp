/**
 * @file iparameters.hpp
 *
 * Defines default initial parameters
 * and reads the input line for changing these parameters
 */
#ifndef INITIALPARAMETERS_HPP
#define INITIALPARAMETERS_HPP

#include <iostream>
#include <string>
#include <getopt.h>
using namespace std;
namespace iparameters {
/*
 Store indices to identify nearest neighbors.
 */
struct InitialParameters  {
    
    InitialParameters();
    
    //Parses the input command line to redefine simulation parameters
    void parse_input(int argc, char** argv);
    
    // prints the basic filename (without extension)
    string base_filename();
    
    
    //Variables in the structure
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
	string order_Q;				// flag defining the system should be initially ordered
	string printC_Q;			// flag defining if the correlation function shall be printed
	string printLastConf_Q;		// flag defining if the last configuration shall be printed
	
};
}
#endif
