#include "iparameters.hpp"

using namespace std;

namespace iparameters {

InitialParameters::InitialParameters(){
	this->Lx = 10;				// in-plane size of the cube lattice
	this->Ly = 10;				// in-plane size of the cube lattice
	this->Lz = 100;				// oo-plane size of the cube lattice
	this->N = Lx * Ly * Lz;		// Number of spins
	this->rs = 0;				// random seed

	this->T = 1.5; 				// Temperature
	this->Hx = 0.0;				// Field strength along z
	this->Hy = 0.0;				// Field strength along z
	this->Hz = 0.0;				// Field strength along z
	
	this->Jcc = 1.0;			// in-plane ferromagnetic coupling (antiferromagnetic if Jxy < 0)
	this->Jab = 1.0;			// oo-plane ferromagnetic coupling (antiferromagnetic if Joo < 0)
	
	this->K = 0.0;				// c-axis anisotropy
	this->L = 0.0;				// in-plane anisotropy
	
	this->eps = 1000.0;						

	this->passes = 1000;		// Default # of passes
	this->printC_Q = "False";					// flag defining if the correlation function shall be printed
	this->printLastConf_Q = "False";			// flag defining if the last configuration shall be printed

	
};

}