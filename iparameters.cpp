#include "iparameters.hpp"
using namespace std;
namespace iparameters {
//constructor
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
	this->order_Q = "False";					// flag defining the system should be initially ordered
	this->printC_Q = "False";					// flag defining if the correlation function shall be printed
	this->printLastConf_Q = "False";			// flag defining if the last configuration shall be printed
	
}
//Parses the input command line to redefine simulation parameters
void InitialParameters::parse_input(int argc, char** argv)
{
	int c;
	while (1)
    {
    	static struct option long_options[] =
        {
        	/* These options don’t set a flag. We distinguish them by their indices. */
            {"T",   required_argument, NULL, 't'},
            {"Hx",  required_argument, NULL, 'x'},
            {"Hy",  required_argument, NULL, 'y'},
            {"Hz",  required_argument, NULL, 'z'},
            {"Jab",  required_argument, NULL, 'j'},
            {"Jcc",  required_argument, NULL, 'm'},
            {"passes",  required_argument, NULL, 'n'},
            {"order",  required_argument, NULL, 'o'},
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
					cout << "error" << endl;		break;
			case 't':
					this->T = stod(optarg);			break;
			case 'x':
					this->Hx = stod(optarg); 		break;
			case 'y':
					this->Hy = stod(optarg); 		break;
			case 'z':
					this->Hz = stod(optarg); 		break;
            case 'j':
                    this->Jab = stod(optarg); 		break;
            case 'm':
                    this->Jcc = stod(optarg); 		break;
            case 'n':
                    this->passes = stoi(optarg); 	break;
            case 'e':
                    this->eps = stod(optarg); 		break;
            case 's':
                    if (soptarg.compare("True") == 0)
                    {
                    this->printLastConf_Q = "True";
                    }
                    else
                    {
                    this->printLastConf_Q = "False";
                    }								break;
            case 'o':
                    if (soptarg.compare("True") == 0)
                    {
                    this->order_Q = "True";
                    }
                    else
                    {
                    this->order_Q = "False";
                    }								break;

            case 'r':
                    this->rs = stoi(optarg);	 	break;
            case '?':
                    /* getopt_long already printed an error message. */ break;
        default:
        abort ();
        }
    }
}
// prints the basic filename (without extension)
string InitialParameters::base_filename(){
    string fn = "Y2-12-7_T=" + to_string(this->T).substr(0,6) + \
    			"_Hz=" + to_string(this->Hz).substr(0,6) + \
    			"_passes=" + to_string(this->passes) + \
    			"_Jab=" + to_string(this->Jab).substr(0,4) + \
    			"_Jcc=" + to_string(this->Jcc).substr(0,4) + \
    			"_eps=" + to_string(this->eps).substr(0,4) + \
    			"_order=" + this->order_Q + \
    			"_rs=" + to_string(this->rs);
	return fn;
}

}