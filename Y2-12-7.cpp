#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <random>
#include <getopt.h>
#include "lattice.hpp"
#include "observables.hpp"
#include "iparameters.hpp"
#include "hamiltonian.hpp"
#include "runningaverage.hpp"
using namespace std;

//Define global parameters

const double pi = M_PI;
double new_acceptance_rate = 0.0;

//Few useful functions to start with:
// random real between 0 and 1
double randomreal() {return (double) rand() / RAND_MAX;}

// random number between -1 and 1
double randomrealpm() {return ((double) 2.0 * randomreal() - 1.0);}

// random integer which assumes either -1 or 1
double randomintpm(){ return (double) 2.0 * round(randomreal()) - 1.0;}

// get vector index from 3d array (i,j,k) \elem [0,Lx-1] * [0,Ly-1] * [0,Lz-1]  ->  l \elem [0, Lx * Ly * Lz - 1] ...
int index(const iparameters::InitialParameters& ip, int i, int j, int k) {return i + j * ip.Lx + k * ip.Lx * ip.Ly;}
// ... and back
int cx(const iparameters::InitialParameters& ip, int l) {return l % ip.Lx;}
int cy(const iparameters::InitialParameters& ip, int l) {return (int) ((l - cx(ip, l)) / ip.Lx) % ip.Ly;}
int cz(const iparameters::InitialParameters& ip, int l) {return (int) ((l - cx(ip, l)) / ip.Lx - cy(ip, l)) / ip.Ly;}

// Initialize a random spin configuration, calculate its energy & magnetization
void initialization(const iparameters::InitialParameters& ip, lattice::EISH::NearestNeighbors& nn, observables::EISH::SystemState& s)
{
	// set the random number generator seed
	srand(ip.rs);
    
    // Generate an initial spin configuration
	// Evaluate for each site the neighboring site indices
	// Define sweep order (initially ordered 0, 1, 2, ...)
	for (int i = 0; i < ip.N; i++)
	{
	    double S;
	    double D;
	    
	    S = randomintpm();
	    D = 0.05 * randomrealpm();
	    s.SA[i] = S;
	    s.DA[i] = D;
	    
	    S = randomintpm();
	    D = 0.05 * randomrealpm();
	    s.SB[i] = S;
	    s.DB[i] = D;
	    	    
        nn.AB1[i] = index(ip, cx(ip, i), cy(ip, i), cz(ip, i));
        nn.AB2[i] = index(ip, (cx(ip, i) + ip.Lx - 1) % ip.Lx, cy(ip, i), cz(ip, i));
        nn.AB3[i] = index(ip, (cx(ip, i) + ip.Lx - 1 + cy(ip, i) % 2) % ip.Lx, (cy(ip, i) + ip.Ly - 1) % ip.Ly, cz(ip, i));
        nn.AB4[i] = index(ip, cx(ip, i), cy(ip, i), (cz(ip, i) + 1) % ip.Lz);
        nn.AB5[i] = index(ip, (cx(ip, i) + ip.Lx - 1) % ip.Lx, cy(ip, i), (cz(ip, i) + 1) % ip.Lz);
        nn.AB6[i] = index(ip, (cx(ip, i) + ip.Lx - 1 + cy(ip, i) % 2) % ip.Lx, (cy(ip, i) + ip.Ly - 1) % ip.Ly, (cz(ip, i) + 1) % ip.Lz);
        nn.AA1[i] = index(ip, cx(ip, i), cy(ip, i), (cz(ip, i) + ip.Lz - 1) % ip.Lz);
        nn.AA2[i] = index(ip, cx(ip, i), cy(ip, i), (cz(ip, i) + 1) % ip.Lz);
        
        nn.BA1[i] = index(ip, (cx(ip, i) + cy(ip, i) % 2) % ip.Lx, (cy(ip, i) + 1) % ip.Ly, (cz(ip, i) + ip.Lz - 1) % ip.Lz);
        nn.BA2[i] = index(ip, cx(ip, i), cy(ip, i), (cz(ip, i) + ip.Lz - 1) % ip.Lz);
        nn.BA3[i] = index(ip, (cx(ip, i) + 1) % ip.Lx, cy(ip, i), (cz(ip, i) + ip.Lz - 1) % ip.Lz);
        nn.BA4[i] = index(ip, (cx(ip, i) + cy(ip, i) % 2) % ip.Lx, (cy(ip, i) + 1) % ip.Ly, cz(ip, i));
        nn.BA5[i] = index(ip, cx(ip, i), cy(ip, i), cz(ip, i));
        nn.BA6[i] = index(ip, (cx(ip, i) + 1) % ip.Lx, cy(ip, i), cz(ip, i));
        nn.BB1[i] = index(ip, cx(ip, i), cy(ip, i), (cz(ip, i) + ip.Lz - 1) % ip.Lz);
        nn.BB2[i] = index(ip, cx(ip, i), cy(ip, i), (cz(ip, i) + 1) % ip.Lz);
        nn.sweeporder[i] = i;
	    
	}

	// Evaluate the Energy and Magnetization of the initial state
	for (int i = 0; i < ip.N; i++)
	{
	    double SA = s.SA[i];
	    double SB = s.SB[i];
	    
	    double DA = s.DA[i];
	    double DB = s.DB[i];
	    
	    s.sysE += hamiltonian::EISH::nnA(ip, SA, DA, i, nn, s)/ 2 + hamiltonian::EISH::onsite(ip, SA, DA, i);
	    s.sysE += hamiltonian::EISH::nnB(ip, SB, DB, i, nn, s)/ 2 + hamiltonian::EISH::onsite(ip, SB, DB, i);

	    s.sysM += SA + SB;
	    s.sysMA += SA;
	    s.sysMB += SB;
	    s.sysDA += DA;
	    s.sysDB += DB;
	    
	    if (cz(ip, i) % 2 == 0)
	    {
	    	s.sysSMA += SA;
	    	s.sysSMB += SB;
	    	s.sysSDA += DA;
	    	s.sysSDB += DB;
	    }
	    else
	    {
	    	s.sysSMA -= SA;
	    	s.sysSMB -= SB;
	    	s.sysSDA -= DA;
	    	s.sysSDB -= DB;
	    }
	    
	    s.sysMAMB += SA * SB;
	    s.sysDADB += DA * DB;
	    s.sysMADA += SA * DA;
	    s.sysMADB += SA * DB;
	    s.sysMBDA += SB * DA;
	    s.sysMBDB += SB * DB;
	    s.sysMAMBDA += SA * SB * DA;
	    s.sysMAMBDB += SA * SB * DB;
	    s.sysMADADB += SA * DA * DB;
	    s.sysMBDADB += SB * DA * DB;
	    s.sysMAMBDADB += SA * SB * DA * DB;
	    
	}
	
	// All quantities are per spin
	s.div_sys_vars(ip.N);
	
        cout << "Energy before ordering : " << s.sysE << endl;
}

// Initialize an ordered spin configuration, calculate its energy & magnetization
void order(const iparameters::InitialParameters& ip, const lattice::EISH::NearestNeighbors& nn, observables::EISH::SystemState& s)
{
    // Generate an initial spin configuration
        // Evaluate for each site the neighboring site indices
        // Define sweep order (initially ordered 0, 1, 2, ...)
        for (int i = 0; i < ip.N; i++)
        {
            if (cz(ip, i) % 2 == 0)
            {
                s.SA[i] = 1.0;
                s.SB[i] = 1.0;
                s.DA[i] = 0.01;
                s.DB[i] = -0.01;
            }
            else
            {
                s.SA[i] = -1.0;
                s.SB[i] = -1.0;
                s.DA[i] = 0.01;
                s.DB[i] = -0.01;
            }
        }


		s.reset_sys_vars();

        // Evaluate the Energy and Magnetization of the ordered state
        for (int i = 0; i < ip.N; i++)
        {
            double SA = s.SA[i];
            double SB = s.SB[i];
            double DA = s.DA[i];
            double DB = s.DB[i];
            
			s.sysE += hamiltonian::EISH::nnA(ip, SA, DA, i, nn, s)/ 2 + hamiltonian::EISH::onsite(ip, SA, DA, i);
			s.sysE += hamiltonian::EISH::nnB(ip, SB, DB, i, nn, s)/ 2 + hamiltonian::EISH::onsite(ip, SB, DB, i);

			s.sysM += SA + SB;
			s.sysMA += SA;
			s.sysMB += SB;
			s.sysDA += DA;
			s.sysDB += DB;
			
			if (cz(ip, i) % 2 == 0)
			{
				s.sysSMA += SA;
				s.sysSMB += SB;
				s.sysSDA += DA;
				s.sysSDB += DB;
			}
			else
			{
				s.sysSMA -= SA;
				s.sysSMB -= SB;
				s.sysSDA -= DA;
				s.sysSDB -= DB;
			}
			
			s.sysMAMB += SA * SB;
			s.sysDADB += DA * DB;
			s.sysMADA += SA * DA;
			s.sysMADB += SA * DB;
			s.sysMBDA += SB * DA;
			s.sysMBDB += SB * DB;
			s.sysMAMBDA += SA * SB * DA;
			s.sysMAMBDB += SA * SB * DB;
			s.sysMADADB += SA * DA * DB;
			s.sysMBDADB += SB * DA * DB;
			s.sysMAMBDADB += SA * SB * DA * DB;
        }
		
		s.div_sys_vars(ip.N);
		
        cout << "Energy after ordering : " << s.sysE << endl;
}

// Create output file and write header line
void fileprint_init(const string& basefn, const string& printline)
{
	ofstream outfile;
	string fn = basefn + ".tsv";
	
	outfile.open(fn);
	outfile << printline + "\n";
	outfile.close();
}


// Export data to file
void fileprint_appendline(const string& basefn, const string& printline)
{
	ofstream outfile;
	string fn = basefn + ".tsv";
	
	outfile.open(fn, ios_base::app);
	outfile << printline + "\n";
	outfile.close();
}


// All necessary functions to generate one pass:
// Generate a random sweep order to visit the spins
void generate_new_sweep_order(const iparameters::InitialParameters& ip, lattice::EISH::NearestNeighbors& nn)
{
	// Shuffle the sweep order
	for (int i = 0; i < ip.N; i++)
	{
	    int mem = nn.sweeporder[i];
	    int newi = random() % ip.N;
	    
	    nn.sweeporder[i] = nn.sweeporder[newi];
	    nn.sweeporder[newi] = mem;
	}
}

// define one spin flip (Metropolis-Hastings step)
// This step also accounts for possible phase space reduction
void oneflipA(const iparameters::InitialParameters& ip, int i, const lattice::EISH::NearestNeighbors& nn, observables::EISH::SystemState& s)
{
	double oldSA = s.SA[i];
	double oldDA = s.DA[i];
	double SB = s.SB[i];
	double DB = s.DB[i];
	
	double newSA = randomintpm();

        // evaluate the maximum absolute displacement of the previous pass
        // and use it as a reference for the next one
        // double newDA = 1.01 * old_max_abs_displacement * randomrealpm();

        double newDA = oldDA + 0.001 * randomrealpm();

	double delta_E = hamiltonian::EISH::siteA(ip, newSA, newDA, i, nn, s) - hamiltonian::EISH::siteA(ip, oldSA, oldDA, i, nn, s);
	
	if ((delta_E < 0) || (randomreal() < exp(- delta_E / ip.T)))
	{
		new_acceptance_rate +=1;
		
        s.SA[i] = newSA;
        s.DA[i] = newDA;
        
        s.sysE += delta_E / ip.N;
        s.sysM += (newSA - oldSA) / ip.N;
        s.sysMA += (newSA - oldSA) / ip.N;
        s.sysDA += (newDA - oldDA) / ip.N;
        
        if (cz(ip, i) % 2 == 0){
        	s.sysSMA += (newSA - oldSA) / ip.N;
        	s.sysSDA += (newDA - oldDA) / ip.N;
        } else {
        	s.sysSMA -= (newSA - oldSA) / ip.N;
        	s.sysSDA -= (newDA - oldDA) / ip.N;
        }
        
        s.sysMAMB += (newSA - oldSA) * SB / ip.N;
		s.sysDADB += (newDA - oldDA) * DB / ip.N;
		s.sysMADA += (newSA * newDA - oldSA * oldDA) / ip.N;
		s.sysMADB += (newSA - oldSA) * DB / ip.N;
		s.sysMBDA += (newDA - oldDA) * SB / ip.N;
		s.sysMAMBDA += (newSA * newDA - oldSA * oldDA) * SB / ip.N;
		s.sysMAMBDB += (newSA - oldSA) * SB * DB / ip.N;
		s.sysMADADB += (newSA * newDA - oldSA * oldDA) * DB / ip.N;
		s.sysMBDADB += (newDA - oldDA) * SB * DB / ip.N;
		s.sysMAMBDADB += (newSA * newDA - oldSA * oldDA) * SB * DB / ip.N;
        
        
    }
}
void oneflipB(const iparameters::InitialParameters& ip, int i, const lattice::EISH::NearestNeighbors& nn, observables::EISH::SystemState& s)
{
	double SA = s.SA[i];
	double DA = s.DA[i];
	double oldSB = s.SB[i];
	double oldDB = s.DB[i];
	
	double newSB = randomintpm();

        //double newDB = 1.01 * old_max_abs_displacement * randomrealpm();

        double newDB = oldDB + 0.001 * randomrealpm();

        double delta_E = hamiltonian::EISH::siteB(ip, newSB, newDB, i, nn, s) - hamiltonian::EISH::siteB(ip, oldSB, oldDB, i, nn, s);
	
	if ((delta_E < 0) || (randomreal() < exp(- delta_E / ip.T)))
	{
		new_acceptance_rate +=1;
		
        s.SB[i] = newSB;
        s.DB[i] = newDB;
        
        s.sysE += delta_E / ip.N;
        s.sysM += (newSB - oldSB) / ip.N;
        s.sysMB += (newSB - oldSB) / ip.N;
        s.sysDB += (newDB - oldDB) / ip.N;
        
        if (cz(ip, i) % 2 == 0){
        	s.sysSMB += (newSB - oldSB) / ip.N;
        	s.sysSDB += (newDB - oldDB) / ip.N;
        } else {
        	s.sysSMB -= (newSB - oldSB) / ip.N;
        	s.sysSDB -= (newDB - oldDB) / ip.N;
        }
        
        s.sysMAMB += SA * (newSB - oldSB) / ip.N;
		s.sysDADB += DA * (newDB - oldDB) / ip.N;
		s.sysMADB += SA * (newDB - oldDB) / ip.N;
		s.sysMBDA += (newSB - oldSB) * DA / ip.N;
		s.sysMBDB += (newSB * newDB - oldSB * oldDB) / ip.N;
		s.sysMAMBDA += SA * (newSB - oldSB) * DA / ip.N;
		s.sysMAMBDB += SA * (newSB * newDB - oldSB * oldDB) / ip.N;
		s.sysMADADB += SA * DA * (newDB - oldDB) / ip.N;
		s.sysMBDADB += DA  * (newSB * newDB - oldSB * oldDB) / ip.N;
		s.sysMAMBDADB += SA * DA * (newSB * newDB - oldSB * oldDB) / ip.N;
        
    } 
}
void oneflip(const iparameters::InitialParameters& ip, int i, const lattice::EISH::NearestNeighbors& nn, observables::EISH::SystemState& s)
{
	oneflipA(ip, i, nn, s);
	oneflipB(ip, ip.N - 1 - i, nn, s);
}

// Defines a full pass of spin flips
// This function evaluates the new acceptance rate and adjusts  the opening angle if necessary
void onepass(const iparameters::InitialParameters& ip, lattice::EISH::NearestNeighbors& nn, observables::EISH::SystemState& s)
{
	generate_new_sweep_order(ip, nn);
	
	new_acceptance_rate = 0.0;

	for (int i = 0; i < ip.N; i++)
	{
        oneflip(ip, nn.sweeporder[i],nn, s);
    }
    
    new_acceptance_rate = new_acceptance_rate / 2 / ip.N;

    s.pass += 1;
}

void compute(iparameters::InitialParameters& ip, lattice::EISH::NearestNeighbors& nn, observables::EISH::SystemState& s, runningaverage::EISH::RunningAverage& averages)
{
	// half of the passes are added as prepasses
	int prepasses = (int) ip.passes / 2;
	// passes are divided into subpasses
	int subpasses = (int) ip.passes/10;
	if(subpasses > 100000){subpasses = 100000;};
	int superpasses = (int) ip.passes / subpasses;
	
	for (int i = 0; i < prepasses; i++)
	{
		onepass(ip, nn, s);
	}
		
	fileprint_init(ip.base_filename(), averages.header_output_string());
	for (int i = 0; i < superpasses; i++)
	{	
		averages.reset();
		
		for (int j = 0; j < subpasses; j++)
		{
			onepass(ip, nn, s);
			averages.increment(s.sysE, s.sysM, s.sysMA, s.sysMB, s.sysSMA, s.sysSMB, s.sysDA, s.sysDB, s.sysSDA, s.sysSDB, s.sysMAMB, s.sysDADB, s.sysMADA, s.sysMADB, s.sysMBDA, s.sysMBDB, s.sysMAMBDA, s.sysMAMBDB, s.sysMADADB, s.sysMBDADB, s.sysMAMBDADB);
		}
		
		averages.div(subpasses);
		
		cout << "Current acceptance rate = " << new_acceptance_rate << endl;

		fileprint_appendline(ip.base_filename(), to_string(s.pass) + "\t" + to_string(ip.T) + "\t" + to_string(ip.Hz) + "\t" + averages.output_string());
		
	}
}
int main(int argc, char** argv)
//-------
{
    iparameters::InitialParameters ip;
    
    ip.parse_input(argc, argv);
    
    cout << "Starting the simulation with the parameters" <<  endl;
    cout << "T = " << ip.T << "   Hz = " << ip.Hz << "   Jab = " << ip.Jab << "   Jcc = " << ip.Jcc << endl;
	time_t t0 = time(nullptr);
    lattice::EISH::NearestNeighbors nn(ip.N);
	observables::EISH::SystemState s(ip.N);
	runningaverage::EISH::RunningAverage averages;
	
	initialization(ip,nn,s);
    if (ip.order_Q == "True") {order(ip, nn, s);}
	compute(ip, nn, s, averages);
	
	time_t t = time(nullptr);
	double dt = 1.0 * (t - t0);
	cout << "Simulation time == " << dt << endl;
	return 0;
}

