#include "lattice.hpp"
using namespace std;
namespace lattice {
namespace EISH {
	NearestNeighbors::NearestNeighbors(int N){
	    this->AB1 = vector<int>(N);
	    this->AB2 = vector<int>(N);
	    this->AB3 = vector<int>(N);
	    this->AB4 = vector<int>(N);
	    this->AB5 = vector<int>(N);
	    this->AB6 = vector<int>(N);
	    this->AA1 = vector<int>(N);
	    this->AA2 = vector<int>(N);
	    this->BA1 = vector<int>(N);
	    this->BA2 = vector<int>(N);
	    this->BA3 = vector<int>(N);
	    this->BA4 = vector<int>(N);
	    this->BA5 = vector<int>(N);
	    this->BA6 = vector<int>(N);
	    this->BB1 = vector<int>(N);
	    this->BB2 = vector<int>(N);
	    this->sweeporder = vector<int>(N);
	}
}    
}
