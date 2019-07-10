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

using namespace std;


//Define global parameters


const double pi = M_PI;
double new_acceptance_rate = 0;


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
void initialization(const iparameters::InitialParameters& ip, lattice::NearestNeighbors& nn, observables::SystemState& s)
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
	    
	    double nDA = s.DA[i];
	    double nDB = s.DB[i];
	        	    
	    s.E += hamiltonian::EISH::nnA(ip, SA, nDA, i, nn, s)/ 2 + hamiltonian::EISH::onsite(ip, SA, nDA, i);
	    
	    s.E += hamiltonian::EISH::nnB(ip, SB, nDB, i, nn, s)/ 2 + hamiltonian::EISH::onsite(ip, SB, nDB, i);
	    
	    s.MA += SA;
	    s.MB += SB;
	    
	    if (cz(ip, i) % 2 == 0)
	    {
	    	s.sMA += SA;
	    	s.sMB += SB;
	    }
	    else
	    {
	    	s.sMA -= SA;
	    	s.sMB -= SB;
	    }
	    
	    s.Lnem += SA * SB;
	    s.Ldel += nDA - nDB;
	    
	}
	
	// All quantities are per spin
	
	s.E /= ip.N;
	s.MA /= ip.N;
	s.MB /= ip.N;
	s.sMA /= ip.N;
	s.sMB /= ip.N;
	s.Lnem /= ip.N;
	s.Ldel /= ip.N;
	
        cout << "Energy before ordering : " << s.E << endl;
}


// Initialize a random spin configuration, calculate its energy & magnetization
void order(const iparameters::InitialParameters& ip, const lattice::NearestNeighbors& nn, observables::SystemState& s)
{
    cout << "Energy before ordering : " << s.E << endl;
    cout << "Distortion before ordering : " << s.Ldel << endl;

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


        s.E = 0.0;
        s.MA = 0.0;
        s.MB = 0.0;
        s.sMA = 0.0;
        s.sMB = 0.0;
        s.Lnem = 0.0;
        s.Ldel = 0.0;


        // Evaluate the Energy and Magnetization of the ordered state
        for (int i = 0; i < ip.N; i++)
        {
            double SA = s.SA[i];
            double SB = s.SB[i];

            double nDA = s.DA[i];
            double nDB = s.DB[i];

            s.E += hamiltonian::EISH::nnA(ip, SA, nDA, i, nn, s)/ 2 + hamiltonian::EISH::onsite(ip, SA, nDA, i);

            s.E += hamiltonian::EISH::nnB(ip, SB, nDB, i, nn, s)/ 2 + hamiltonian::EISH::onsite(ip, SB, nDB, i);

            s.MA += SA;
            s.MB += SB;

            if (cz(ip, i) % 2 == 0)
            {
                s.sMA += SA;
                s.sMB += SB;
            }
            else
            {
                s.sMA -= SA;
                s.sMB -= SB;
            }

            s.Lnem += SA * SB;
            s.Ldel += nDA - nDB;

        }

        // All quantities are per spin

        s.E /= ip.N;
        s.MA /= ip.N;
        s.MB /= ip.N;
        s.sMA /= ip.N;
        s.sMB /= ip.N;
        s.Lnem /= ip.N;
        s.Ldel /= ip.N;

        cout << "Energy after ordering : " << s.E << endl;
        cout << "Distortion after ordering : " << s.Ldel << endl;
}


//All necessary functions for printing out useful output
//Generate filename for the main data file
string filename(const iparameters::InitialParameters& ip)
{
    string fn = "Y2-12-7_T=" + to_string(ip.T).substr(0,6) + "_Hz=" + to_string(ip.Hz).substr(0,6) + "_passes=" + to_string(ip.passes) + "_Jab=" + to_string(ip.Jab).substr(0,4) + "_eps=" + to_string(ip.eps).substr(0,4) + "_Jcc=" + to_string(ip.Jcc).substr(0,4) + "_rs=" + to_string(ip.rs);
	return fn;
}



// Create output file and write header line
void init_printf(const iparameters::InitialParameters& ip)
{
	ofstream outfile;
	string fn = filename(ip) + ".tsv";
	
	outfile.open(fn);
	outfile << "pass \t Temperature \t Hz \t <E> \t <E^2> \t <MA> \t <MA^2> \t <MB> \t <MB^2> \t <sMA> \t <sMA^2> \t <sMB> \t <sMB^2> \t <Lnem> \t <Lnem^2> \t <Ldel> \t <Ldel^2> \n";
	outfile.close();
}

// Export data to file
void printf(const iparameters::InitialParameters& ip, const observables::SystemState& s, double e, double e2, double ma, double ma2, double mb, double mb2, double sma, double sma2, double smb, double smb2, double lnem, double lnem2, double ldel, double ldel2)
{
	ofstream outfile;
	string fn = filename(ip) + ".tsv";
	
	outfile.open(fn, ios_base::app);
	outfile << to_string(s.pass) + "\t" + to_string(ip.T) + "\t" + to_string(ip.Hz) + "\t" + to_string(e) + "\t" + to_string(e2) + "\t" + to_string(ma) + "\t" + to_string(ma2) + "\t" + to_string(mb) + "\t" + to_string(mb2) + "\t" + to_string(sma) + "\t" + to_string(sma2) + "\t" + to_string(smb) + "\t" + to_string(smb2) + "\t" + to_string(lnem) + "\t" + to_string(lnem2) + "\t" + to_string(ldel) + "\t" + to_string(ldel2) + "\n";
	outfile.close();
}



//Generate filename for the last configuration data file
string filenameLSC(const iparameters::InitialParameters& ip)
{
        string fn = filename(ip) + "_spinstate" + ".tsv";
        return fn;
}

string filenameLDC(const iparameters::InitialParameters& ip)
{
        string fn = filename(ip) + "_displacement" + ".tsv";
        return fn;
}

// Create output file for last configuration
void init_printLSC(const iparameters::InitialParameters& ip)
{
        ofstream outfile;
        string fn = filenameLSC(ip);

        outfile.open(fn);
        outfile << "Spinstate of a " << ip.Lx << "x" << ip.Ly << "x" << ip.Lz << " Ising spin system  szA0, szA1, ...  dann szB0, szB1, ... \n";
        outfile.close();
}

void init_printLDC(const iparameters::InitialParameters& ip)
{
        ofstream outfile;
        string fn = filenameLDC(ip);

        outfile.open(fn);
        outfile << "Displacement state of a " << ip.Lx << "x" << ip.Ly << "x" << ip.Lz << " Ising spin system  dz0, dz1, ... dann dzB0, dzB1, ... \n";
        outfile.close();
}

void printLSC(const iparameters::InitialParameters& ip, const observables::SystemState& s)
{
        ofstream outfile;
        string fn = filenameLSC(ip);

        outfile.open(fn,ios_base::app);
        for (int i = 0; i < ip.N; i++)
        {
                outfile << to_string(s.SA[i]) << "\n";
        }

        outfile << "\n";

        for (int i = 0; i < ip.N; i++)
        {
                outfile << to_string(s.SB[i]) << "\n";
        }
        outfile.close();
}

void printLDC(const iparameters::InitialParameters& ip, const observables::SystemState& s)
{
        ofstream outfile;
        string fn = filenameLDC(ip);

        outfile.open(fn,ios_base::app);
        for (int i = 0; i < ip.N; i++)
        {
                outfile << to_string(s.DA[i]) << "\n";
        }

        outfile << "\n";

        for (int i = 0; i < ip.N; i++)
        {
                outfile << to_string(s.DB[i]) << "\n";
        }
        outfile.close();
}





























// All necessary functions to generate one pass:
// Generate a random sweep order to visit the spins
void generate_new_sweep_order(const iparameters::InitialParameters& ip, lattice::NearestNeighbors& nn)
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
void oneflipA(const iparameters::InitialParameters& ip, int i, const lattice::NearestNeighbors& nn, observables::SystemState& s)
{
	double SA = s.SA[i];
	double oDA = s.DA[i];
	
	double newSA = randomintpm();
	//double newDA = oDA * ( 1.0 + 3.0 * randomrealpm()) / 2.0;
	double newDA = oDA + 0.01 * randomrealpm();
	double delta_E = hamiltonian::EISH::siteA(ip, newSA, newDA, i, nn, s) - hamiltonian::EISH::siteA(ip, SA, oDA, i, nn, s);
	
	if ((delta_E < 0) || (randomreal() < exp(- delta_E / ip.T)))
	{
		new_acceptance_rate +=1;
		
        s.SA[i] = newSA;
        s.DA[i] = newDA;
        
        s.E += delta_E / ip.N;
        s.MA += (newSA - SA) / ip.N;
        
        if (cz(ip, i) % 2 == 0) {s.sMA += (newSA - SA) / ip.N;}
	    else {s.sMA -= (newSA - SA) / ip.N;}
        
        s.Lnem += (newSA - SA) * s.SB[i] / ip.N;
        s.Ldel += (newDA - oDA) / ip.N;
		
    } 
}

void oneflipB(const iparameters::InitialParameters& ip, int i, const lattice::NearestNeighbors& nn, observables::SystemState& s)
{
	double SB = s.SB[i];
	double oDB = s.DB[i];
	
	double newSB = randomintpm();
	//double newDB = oDB * ( 1.0 + 3.0 * randomrealpm()) / 2.0;
	double newDB = oDB + 0.01 * randomrealpm();
	double delta_E = hamiltonian::EISH::siteB(ip, newSB, newDB, i, nn, s) - hamiltonian::EISH::siteB(ip, SB, oDB, i, nn, s);
	
	if ((delta_E < 0) || (randomreal() < exp(- delta_E / ip.T)))
	{
		new_acceptance_rate +=1;
		
        s.SB[i] = newSB;
        s.DB[i] = newDB;
        
        s.E += delta_E / ip.N;
        s.MB += (newSB - SB) / ip.N;
        
        if (cz(ip, i) % 2 == 0) {s.sMB += (newSB - SB) / ip.N;}
	    else {s.sMB -= (newSB - SB) / ip.N;}
        
        s.Lnem += (newSB - SB) * s.SA[i] / ip.N;
		s.Ldel -= (newDB - oDB) / ip.N;
    } 
}

void oneflip(const iparameters::InitialParameters& ip, int i, const lattice::NearestNeighbors& nn, observables::SystemState& s)
{
	oneflipA(ip, i, nn, s);
	oneflipB(ip, ip.N - 1 - i, nn, s);
}



// Defines a full pass of spin flips
// This function evaluates the new acceptance rate and adjusts  the opening angle if necessary
void onepass(const iparameters::InitialParameters& ip, lattice::NearestNeighbors& nn, observables::SystemState& s)
{
	generate_new_sweep_order(ip, nn);
	
	new_acceptance_rate = 0.0;
	
	for (int i = 0; i < ip.N; i++)
	{
        oneflip(ip, nn.sweeporder[i],nn, s);
    }
    
    new_acceptance_rate = new_acceptance_rate / 2 / ip.N;
}


void compute(const iparameters::InitialParameters& ip, lattice::NearestNeighbors& nn, observables::SystemState& s)
{
	// roughly half of the passes should be added as prepasses
	int prepasses = (int) ip.passes / 2;
	if (prepasses > 500000)
		prepasses = 500000;
	
	int subpasses = 100000;
	if (ip.passes <= subpasses){subpasses = (int) ip.passes/10;}
	
	int superpasses = (int) ip.passes / subpasses;
	
	
	s.pass = 0;
	for (int i = 0; i < prepasses; i++)
	{
		onepass(ip, nn,s);
		s.pass += 1;
		
	}
	
	init_printf(ip);

	for (int i = 0; i < superpasses; i++)
	{	
		double e = 0.0;
		double e2 = 0.0;
		double ma = 0.0;
		double mb = 0.0;
		double ma2 = 0.0;
		double mb2 = 0.0;
		double sma = 0.0;
		double smb = 0.0;
		double sma2 = 0.0;
		double smb2 = 0.0;
		double lnem = 0.0;
		double ldel = 0.0;
		double lnem2 = 0.0;
		double ldel2 = 0.0;
		
		for (int j = 0; j < subpasses; j++)
		{
			onepass(ip, nn,s);
			e += s.E;
			e2 += s.E * s.E;

			ma += s.MA;
			mb += s.MB;
			
			ma2 += s.MA * s.MA;
			mb2 += s.MB * s.MB;
			
			sma += s.sMA;
			smb += s.sMB;
			
			sma2 += s.sMA * s.sMA;
			smb2 += s.sMB * s.sMB;
			
			
			lnem += s.Lnem;
			ldel += s.Ldel;
			
			lnem2 += s.Lnem * s.Lnem;
			ldel2 += s.Ldel * s.Ldel;
			
			s.pass += 1;
		}
		
		e /= subpasses;
		e2 /= subpasses;
	
		ma /= subpasses;
		mb /= subpasses;
		
		ma2 /= subpasses;
		mb2 /= subpasses;
			
		sma /= subpasses;
		smb /= subpasses;
		
		sma2 /= subpasses;
		smb2 /= subpasses;
		
		lnem /= subpasses;
		ldel /= subpasses;
		
		lnem2 /= subpasses;
		ldel2 /= subpasses;
		
		cout << "Current acceptance rate = " << new_acceptance_rate << endl;
		
		printf(ip, s, e, e2, ma, ma2, mb, mb2, sma, sma2, smb, smb2, lnem, lnem2, ldel, ldel2);
		
	}
	
        if (ip.printLastConf_Q.compare("True") == 0)
	{
                init_printLSC(ip);
                printLSC(ip, s);
                init_printLDC(ip);
                printLDC(ip, s);
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

    lattice::NearestNeighbors nn(ip.N);
	observables::SystemState s(ip.N);
	
	initialization(ip,nn,s);

    order(ip, nn, s);

	compute(ip, nn, s);
	
	time_t t = time(nullptr);
	double dt = 1.0 * (t - t0);
	cout << "Simulation time == " << dt << endl;
	return 0;
}


