#include "observables.hpp"

namespace observables{

using namespace std;

SystemState::SystemState(int N){
    //initialize structure elements
    this->E = 0.0;							// Energy per spin

    this->MA = 0.0;							// Magnetization of sublattice A
    this->MB = 0.0;							// Magnetization of sublattice B

    this->sMA = 0.0;						// Staggered Magnetization A
    this->sMB = 0.0;						// Staggered Magnetization B
        
    this->Lnem = 0.0;						// nematic order SA * SB
    this->Ldel = 0.0;						// DA - DB

    this->pass = 0;							// current pass
	
	this->SA = vector<int>(N);
    this->SB = vector<int>(N);
    this->DA = vector<float>(N);
    this->DB = vector<float>(N);
};
	
void SystemState::print_E(){
	cout << this->E << endl;
};

}