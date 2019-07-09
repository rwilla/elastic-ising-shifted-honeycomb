#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <random>
#include <getopt.h>
#include "lattice.cpp"


using namespace std;


//Define global parameters
const int xysize = 10;						// in-plane size of the cube lattice
const int oosize = 100;						// oo-plane size of the cube lattice
const int N = xysize * xysize * oosize;		// Number of spins
const double pi = M_PI;						// pi
const double thr_acceptance_rate = 0.5;		// threshold acceptance rate
int rs = 0;									// random seed

double T = 1.5; 							// Temperature
double Hz = 0.0;							// Field strength along z
int passes = 1000;							// Default # of passes
int pass = 0;								// pass number
double Jcc = 1.0;							// in-plane ferromagnetic coupling (antiferromagnetic if Jxy < 0)
double Jab = 1.0;							// oo-plane ferromagnetic coupling (antiferromagnetic if Joo < 0)
double eps = 1000.0;						// in-plane anisotropy
double epsover2 = eps / 2.0;				// in-plane anisotropy

struct SystemState {
        SystemState(int N);

        double E = 0.0;								// Energy per spin

        double MA = 0.0;							// Magnetization of sublattice A
        double MB = 0.0;							// Magnetization of sublattice B

        double sMA = 0.0;							// Staggered Magnetization A
        double sMB = 0.0;							// Staggered Magnetization B

        double Lnem = 0.0;							// nematic order MA * MB
        double Ldel = 0.0;							// DA - DB

        //double Cfunc [(int) xysize / 2];			// Correlation function <S(0)S(r)>

        std::vector<int> SA;							// Spin state S in sublattice A
        std::vector<int> SB;							// Spin state S in sublattice B
        std::vector<int> DA;							// Distortion state in sublattice A
        std::vector<int> DB;							// Distortion state in sublattice B
};

SystemState::SystemState(int N){
        this->SA = std::vector<int>(N);
        this->SB = std::vector<int>(N);
        this->DA = std::vector<int>(N);
        this->DB = std::vector<int>(N);
}


double displacement_amplitude = .25;		// maximal displacement amplitude
double acceptance_rate = 1.0;				// acceptance rate of the finished pass
double new_acceptance_rate = 1.0;			// new acceptance rate of the running pass
double av_acceptance_rate = 0.0;			// average acceptance rate over all subpasses of one superpass


string printC_Q = "False";					// flag defining if the correlation function shall be printed
string printLastConf_Q = "False";			// flag defining if the last configuration shall be printed




//Parses the input command line to redefine simulation parameters
void parse_input(int argc, char** argv)
{
	int c;
	while (1)
    {
    	static struct option long_options[] =
        {
        	/* These options don’t set a flag. We distinguish them by their indices. */
            {"T",   required_argument, NULL, 't'},
            {"Hz",  required_argument, NULL, 'z'},
            {"Jab",  required_argument, NULL, 'j'},
            {"Jcc",  required_argument, NULL, 'm'},
            {"passes",  required_argument, NULL, 'n'},
            {"eps",  required_argument, NULL, 'e'},
            //{"printCfunction",  required_argument, NULL, 'c'},
            {"printSpinConfiguration",  required_argument, NULL, 's'},
            {"rs",  required_argument, NULL, 'r'},
            {NULL, 0, NULL, 0}
        };
    	/* getopt_long stores the option index here. */
    	int option_index = 0;

		c = getopt_long_only (argc, argv, "a:b:c:d:", long_options, NULL);

		/* Detect the end of the options. */
		if (c == -1)
        	break;
        	
        string soptarg = optarg;
        switch (c)
        {
        	case 0:
                    cout << "error" << endl;	break;
            case 't':
                    T = stod(optarg);		break;
            case 'z':
                    Hz = stod(optarg); 		break;
            case 'j':
                    Jab = stod(optarg);	 	break;
            case 'm':
                    Jcc = stod(optarg); 	break;
            case 'n':
                    passes = stoi(optarg); 		break;
            case 'e':
                    {eps = stod(optarg); epsover2 = eps/2.0;} 		break;
            /*
            case 'c':
                    if (soptarg.compare("True") == 0)
                    {
                            printC_Q = "True";
                    } else {
                            printC_Q = "False";
                    }				break;
            */
            case 's':
                    if (soptarg.compare("True") == 0)
                    {
                            printLastConf_Q = "True";
                    } else {
                            printLastConf_Q = "False";
                    }				break;
            case 'r':
                    rs = stoi(optarg);	 	break;
            case '?':
                    /* getopt_long already printed an error message. */ break;
        default:
        abort ();
        }
    }
}







//Few useful functions to start with:
// random real between 0 and 1
double randomreal() {return (double) rand() / RAND_MAX;}

// random number between -1 and 1
double randomrealpm() {return ((double) 2.0 * randomreal() - 1.0);}

// random integer which assumes either -1 or 1
double randomintpm(){ return (double) 2.0 * round(randomreal()) - 1.0;}

// get vector index from 3d array (i,j,k) \elem [0,xysize-1]^2 * [0,oosize-1]  ->  l \elem [0, xysize^2 * oosize - 1] ...
int index(int i, int j, int k) {return i + j * xysize + k * xysize * xysize;}

// ... and back
int cx(int l) {return l % xysize;}
int cy(int l) {return (int) ((l - cx(l)) / xysize) % xysize;}
int cz(int l) {return (int) ((l - cx(l)) / xysize - cy(l)) / xysize;}


// Define 3D Anisotropic Heisenberg model with onsite and nn interactions:
// Onsite energy
double hamiltonian_onsite(double S, double D, int spin_index)
{
    double e = 0.0;
    
    //displacement energy
    e += epsover2 * D * D;

    //field-dependence
    e -= Hz * S;
    
    return e;
}


//Nearest neighbor interaction for a spin on sublattice A
double hamiltonian_nnA(double S, double D, int spin_index, const lattice::NearestNeighbors& nn, SystemState& s)
{
    double e = 0.0;
        
    double nSB;
    double nDB;
            
    nSB = s.SB[nn.AB1[spin_index]];
    nDB = s.DB[nn.AB1[spin_index]];
    
    e += - Jab * (1.0 + D - nDB) * (S * nSB);
    
    nSB = s.SB[nn.AB2[spin_index]];
    nDB = s.DB[nn.AB2[spin_index]];
    
    e += - Jab * (1.0 + D - nDB) * (S * nSB);
    
    nSB = s.SB[nn.AB3[spin_index]];
    nDB = s.DB[nn.AB3[spin_index]];
    
    e += - Jab * (1.0 + D - nDB) * (S * nSB);
    
    nSB = s.SB[nn.AB4[spin_index]];
    nDB = s.DB[nn.AB4[spin_index]];
    
    e += - Jab * (1.0 - D + nDB) * (S * nSB);
    
    nSB = s.SB[nn.AB5[spin_index]];
    nDB = s.DB[nn.AB5[spin_index]];
    
    e += - Jab * (1.0 - D + nDB) * (S * nSB);
    
    nSB = s.SB[nn.AB6[spin_index]];
    nDB = s.DB[nn.AB6[spin_index]];
    
    e += - Jab * (1.0 - D + nDB) * (S * nSB);
    
    double nSA;
    double nDA;

    nSA = s.SA[nn.AA1[spin_index]];
    nDA = s.DA[nn.AA1[spin_index]];
    
    e += - Jcc * (1.0 + D - nDA) * (S * nSA);
    
    nSA = s.SA[nn.AA2[spin_index]];
    nDA = s.DA[nn.AA2[spin_index]];
    
    e += - Jcc * (1.0 - D + nDA) * (S * nSA);
    
    return e;
}


//Nearest neighbor interaction for a spin on sublattice B
double hamiltonian_nnB(double S, double D, int spin_index, const lattice::NearestNeighbors& nn, SystemState& s)
{
    double e = 0.0;
        
    double nSA;
    double nDA;
            
    nSA = s.SA[nn.BA1[spin_index]];
    nDA = s.DA[nn.BA1[spin_index]];
    
    e += - Jab * (1.0 + D - nDA) * (S * nSA);
    
    nSA = s.SA[nn.BA2[spin_index]];
    nDA = s.DA[nn.BA2[spin_index]];
    
    e += - Jab * (1.0 + D - nDA) * (S * nSA);
    
    nSA = s.SA[nn.BA3[spin_index]];
    nDA = s.DA[nn.BA3[spin_index]];
    
    e += - Jab * (1.0 + D - nDA) * (S * nSA);
    
    nSA = s.SA[nn.BA4[spin_index]];
    nDA = s.DA[nn.BA4[spin_index]];
    
    e += - Jab * (1.0 - D + nDA) * (S * nSA);
    
    nSA = s.SA[nn.BA5[spin_index]];
    nDA = s.DA[nn.BA5[spin_index]];
    
    e += - Jab * (1.0 - D + nDA) * (S * nSA);
    
    nSA = s.SA[nn.BA6[spin_index]];
    nDA = s.DA[nn.BA6[spin_index]];
    
    e += - Jab * (1.0 - D + nDA) * (S * nSA);
    
    double nSB;
    double nDB;

    nSB = s.SB[nn.BB1[spin_index]];
    nDB = s.DB[nn.BB1[spin_index]];
    
    e += - Jcc * (1.0 + D - nDB) * (S * nSB);
    
    nSB = s.SB[nn.BB2[spin_index]];
    nDB = s.DB[nn.BB2[spin_index]];
    
    e += - Jcc * (1.0 - D + nDB) * (S * nSB);
    
    return e;
}


//Full hamiltonian for a spin on sublattice A
double hamiltonianA(double S, double D, int spin_index, const lattice::NearestNeighbors& nn, SystemState& s)
{
    double e = 0.0;
    e += hamiltonian_onsite(S, D, spin_index);
    e += hamiltonian_nnA(S, D, spin_index, nn, s);
    
    return e;
}

//Full hamiltonian for a spin on sublattice B
double hamiltonianB(double S, double D, int spin_index, const lattice::NearestNeighbors& nn, SystemState& s)
{
    double e = 0.0;
    e += hamiltonian_onsite(S, D, spin_index);
    e += hamiltonian_nnB(S, D, spin_index, nn, s);
    
    return e;
}


// Initialize a random spin configuration, calculate its energy & magnetization
void initialization(lattice::NearestNeighbors& nn, SystemState& s)
{
	// set the random number generator seed
	srand(rs);
    
    // Generate an initial spin configuration
	// Evaluate for each site the neighboring site indices
	// Define sweep order (initially ordered 0, 1, 2, ...)
	for (int i = 0; i < N; i++)
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
	    	    
        nn.AB1[i] = index(cx(i), cy(i), cz(i));
        nn.AB2[i] = index((cx(i) + xysize - 1) % xysize, cy(i), cz(i));
        nn.AB3[i] = index((cx(i) + xysize - 1 + cy(i) % 2) % xysize, (cy(i) + xysize - 1) % xysize, cz(i));
        nn.AB4[i] = index(cx(i), cy(i), (cz(i) + 1) % oosize);
        nn.AB5[i] = index((cx(i) + xysize - 1) % xysize, cy(i), (cz(i) + 1) % oosize);
        nn.AB6[i] = index((cx(i) + xysize - 1 + cy(i) % 2) % xysize, (cy(i) + xysize - 1) % xysize, (cz(i) + 1) % oosize);

		nn.AA1[i] = index(cx(i), cy(i), (cz(i) + oosize - 1) % oosize);
        nn.AA2[i] = index(cx(i), cy(i), (cz(i) + 1) % oosize);
        
        nn.BA1[i] = index((cx(i) + cy(i) % 2) % xysize, (cy(i) + 1) % xysize, (cz(i) + oosize - 1) % oosize);
        nn.BA2[i] = index(cx(i), cy(i), (cz(i) + oosize - 1) % oosize);
        nn.BA3[i] = index((cx(i) + 1) % xysize, cy(i), (cz(i) + oosize - 1) % oosize);
        nn.BA4[i] = index((cx(i) + cy(i) % 2) % xysize, (cy(i) + 1) % xysize, cz(i));
        nn.BA5[i] = index(cx(i), cy(i), cz(i));
        nn.BA6[i] = index((cx(i) + 1) % xysize, cy(i), cz(i));

		nn.BB1[i] = index(cx(i), cy(i), (cz(i) + oosize - 1) % oosize);
        nn.BB2[i] = index(cx(i), cy(i), (cz(i) + 1) % oosize);

        nn.sweeporder[i] = i;
	    
	}


	// Evaluate the Energy and Magnetization of the initial state
	for (int i = 0; i < N; i++)
	{
	    double SA = s.SA[i];
	    double SB = s.SB[i];
	    
	    double nDA = s.DA[i];
	    double nDB = s.DB[i];
	        	    
	    s.E += hamiltonian_nnA(SA, nDA, i, nn, s)/ 2 + hamiltonian_onsite(SA, nDA, i);
	    
	    s.E += hamiltonian_nnB(SB, nDB, i, nn, s)/ 2 + hamiltonian_onsite(SB, nDB, i);
	    
	    s.MA += SA;
	    s.MB += SB;
	    
	    if (cz(i) % 2 == 0)
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
	
	s.E /= N;
	s.MA /= N;
	s.MB /= N;
	s.sMA /= N;
	s.sMB /= N;
	s.Lnem /= N;
	s.Ldel /= N;
	
        cout << "Energy before ordering : " << s.E << endl;
}


// Initialize a random spin configuration, calculate its energy & magnetization
void order(const lattice::NearestNeighbors& nn, SystemState& s)
{
    cout << "Energy before ordering : " << s.E << endl;
    cout << "Distortion before ordering : " << s.Ldel << endl;

    // Generate an initial spin configuration
        // Evaluate for each site the neighboring site indices
        // Define sweep order (initially ordered 0, 1, 2, ...)
        for (int i = 0; i < N; i++)
        {
            if (cz(i) % 2 == 0)
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
        for (int i = 0; i < N; i++)
        {
            double SA = s.SA[i];
            double SB = s.SB[i];

            double nDA = s.DA[i];
            double nDB = s.DB[i];

            s.E += hamiltonian_nnA(SA, nDA, i, nn, s)/ 2 + hamiltonian_onsite(SA, nDA, i);

            s.E += hamiltonian_nnB(SB, nDB, i, nn, s)/ 2 + hamiltonian_onsite(SB, nDB, i);

            s.MA += SA;
            s.MB += SB;

            if (cz(i) % 2 == 0)
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

        s.E /= N;
        s.MA /= N;
        s.MB /= N;
        s.sMA /= N;
        s.sMB /= N;
        s.Lnem /= N;
        s.Ldel /= N;

        cout << "Energy after ordering : " << s.E << endl;
        cout << "Distortion after ordering : " << s.Ldel << endl;
}


//All necessary functions for printing out useful output
//Generate filename for the main data file
string filename()
{
    string fn = "Y2-12-7_T=" + to_string(T).substr(0,6) + "_Hz=" + to_string(Hz).substr(0,6) + "_passes=" + to_string(passes) + "_Jab=" + to_string(Jab).substr(0,4) + "_eps=" + to_string(eps).substr(0,4) + "_Jcc=" + to_string(Jcc).substr(0,4) + "_rs=" + to_string(rs);
	return fn;
}



// Create output file and write header line
void init_printf()
{
	ofstream outfile;
	string fn = filename() + ".tsv";
	
	outfile.open(fn);
	outfile << "pass \t Temperature \t Hz \t <E> \t <E^2> \t <MA> \t <MA^2> \t <MB> \t <MB^2> \t <sMA> \t <sMA^2> \t <sMB> \t <sMB^2> \t <Lnem> \t <Lnem^2> \t <Ldel> \t <Ldel^2> \n";
	outfile.close();
}

// Export data to file
void printf(double e, double e2, double ma, double ma2, double mb, double mb2, double sma, double sma2, double smb, double smb2, double lnem, double lnem2, double ldel, double ldel2)
{
	ofstream outfile;
	string fn = filename() + ".tsv";
	
	outfile.open(fn,ios_base::app);
	outfile << to_string(pass) + "\t" + to_string(T) + "\t" + to_string(Hz) + "\t" + to_string(e) + "\t" + to_string(e2) + "\t" + to_string(ma) + "\t" + to_string(ma2) + "\t" + to_string(mb) + "\t" + to_string(mb2) + "\t" + to_string(sma) + "\t" + to_string(sma2) + "\t" + to_string(smb) + "\t" + to_string(smb2) + "\t" + to_string(lnem) + "\t" + to_string(lnem2) + "\t" + to_string(ldel) + "\t" + to_string(ldel2) + "\n";
	outfile.close();
}



//Generate filename for the last configuration data file
string filenameLSC()
{
        string fn = filename() + "_spinstate" + ".tsv";
        return fn;
}

string filenameLDC()
{
        string fn = filename() + "_displacement" + ".tsv";
        return fn;
}

// Create output file for last configuration
void init_printLSC()
{
        ofstream outfile;
        string fn = filenameLSC();

        outfile.open(fn);
        outfile << "Spinstate of a " << xysize << "x" << xysize << "x" << oosize << " Ising spin system  szA0, szA1, ...  dann szB0, szB1, ... \n";
        outfile.close();
}

void init_printLDC()
{
        ofstream outfile;
        string fn = filenameLDC();

        outfile.open(fn);
        outfile << "Displacement state of a " << xysize << "x" << xysize << "x" << oosize << " Ising spin system  dz0, dz1, ... dann dzB0, dzB1, ... \n";
        outfile.close();
}

void printLSC(const SystemState& s)
{
        ofstream outfile;
        string fn = filenameLSC();

        outfile.open(fn,ios_base::app);
        for (int i = 0; i < N; i++)
        {
                outfile << to_string(s.SA[i]) << "\n";
        }

        outfile << "\n";

        for (int i = 0; i < N; i++)
        {
                outfile << to_string(s.SB[i]) << "\n";
        }
        outfile.close();
}

void printLDC(const SystemState& s)
{
        ofstream outfile;
        string fn = filenameLDC();

        outfile.open(fn,ios_base::app);
        for (int i = 0; i < N; i++)
        {
                outfile << to_string(s.DA[i]) << "\n";
        }

        outfile << "\n";

        for (int i = 0; i < N; i++)
        {
                outfile << to_string(s.DB[i]) << "\n";
        }
        outfile.close();
}





























// All necessary functions to generate one pass:
// Generate a random sweep order to visit the spins
void generate_new_sweep_order(lattice::NearestNeighbors& nn)
{
	// Shuffle the sweep order
	for (int i = 0; i < N; i++)
	{
	    int mem = nn.sweeporder[i];
	    int newi = random() % N;
	    
	    nn.sweeporder[i] = nn.sweeporder[newi];
	    nn.sweeporder[newi] = mem;
	}
}


// define one spin flip (Metropolis-Hastings step)
// This step also accounts for possible phase space reduction
void oneflipA(int i, const lattice::NearestNeighbors& nn, SystemState& s)
{
	double SA = s.SA[i];
	double oDA = s.DA[i];
	
	double newSA = randomintpm();
	//rw- 0.01
	double newDA = oDA + 0.01 * randomrealpm();
		
	double delta_E = hamiltonianA(newSA, newDA, i, nn, s) - hamiltonianA(SA, oDA, i, nn, s);
	
	if ((delta_E < 0) || (randomreal() < exp(- delta_E / T)))
	{
		new_acceptance_rate +=1;
		
        s.SA[i] = newSA;
        s.DA[i] = newDA;
        
        s.E += delta_E / N;
        s.MA += (newSA - SA) / N;
        
        if (cz(i) % 2 == 0) {s.sMA += (newSA - SA) / N;}
	    else {s.sMA -= (newSA - SA) / N;}
        
        s.Lnem += (newSA - SA) * s.SB[i] / N;
        s.Ldel += (newDA - oDA) / N;
		
    } 
}

void oneflipB(int i, const lattice::NearestNeighbors& nn, SystemState& s)
{
	double SB = s.SB[i];
	double oDB = s.DB[i];
	
	double newSB = randomintpm();
	//rw- 0.01
	double newDB = oDB + 0.01 * randomrealpm();
		
	double delta_E = hamiltonianB(newSB, newDB, i, nn, s) - hamiltonianB(SB, oDB, i, nn, s);
	
	if ((delta_E < 0) || (randomreal() < exp(- delta_E / T)))
	{
		new_acceptance_rate +=1;
		
        s.SB[i] = newSB;
        s.DB[i] = newDB;
        
        s.E += delta_E / N;
        s.MB += (newSB - SB) / N;
        
        if (cz(i) % 2 == 0) {s.sMB += (newSB - SB) / N;}
	    else {s.sMB -= (newSB - SB) / N;}
        
        s.Lnem += (newSB - SB) * s.SA[i] / N;
		s.Ldel -= (newDB - oDB) / N;
    } 
}

void oneflip(int i, const lattice::NearestNeighbors& nn, SystemState& s)
{
	oneflipA(i, nn, s);
	oneflipB(N - 1 - i, nn, s);
}



// Defines a full pass of spin flips
// This function evaluates the new acceptance rate and adjusts  the opening angle if necessary
void onepass(lattice::NearestNeighbors& nn, SystemState& s)
{
	generate_new_sweep_order(nn);
	
	new_acceptance_rate = 0.0;
	
	for (int i = 0; i < N; i++)
	{
        oneflip(nn.sweeporder[i],nn, s);
    }
    
    new_acceptance_rate = new_acceptance_rate / 2 / N;
}


void compute(lattice::NearestNeighbors& nn, SystemState& s)
{
	// roughly half of the passes should be added as prepasses
	int prepasses = (int) passes / 2;
	if (prepasses > 500000)
		prepasses = 500000;
	
	int subpasses = 100000;
	if (passes <= subpasses){subpasses = (int) passes/10;}
	
	int superpasses = (int) passes / subpasses;
	
	
	pass = 0;
	for (int i = 0; i < prepasses; i++)
	{
		onepass(nn,s);
		pass += 1;
		
	}
	
	init_printf();

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
			onepass(nn,s);
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
			
			pass += 1;
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
		
		printf(e, e2, ma, ma2, mb, mb2, sma, sma2, smb, smb2, lnem, lnem2, ldel, ldel2);
		
	}
	
        if (printLastConf_Q.compare("True") == 0)
	{
                init_printLSC();
                printLSC(s);
                init_printLDC();
                printLDC(s);
	}

}

int main(int argc, char** argv)
//-------
{
    parse_input(argc, argv);
    
    cout << "Starting the simulation with the parameters" <<  endl;
    cout << "T = " << T << "   Hz = " << Hz << "   Jab = " << Jab << "   Jcc = " << Jcc << endl;
	time_t t0 = time(nullptr);

        lattice::NearestNeighbors nn(N);
        SystemState s(N);
	
	initialization(nn,s);

        order(nn, s);

	compute(nn, s);
	
	time_t t = time(nullptr);
	double dt = 1.0 * (t - t0);
	cout << "Simulation time == " << dt << endl;
	return 0;
}


