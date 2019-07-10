#include <iostream>
#include "observables.hpp"

namespace observables{

SystemState::SystemState(int N){
    //initialize structure elements
    this->E = 0.0;							// Energy per spin

    this->MA = 0.0;							// Magnetization of sublattice A
    this->MB = 0.0;							// Magnetization of sublattice B

    this->sMA = 0.0;						// Staggered Magnetization A
    this->sMB = 0.0;						// Staggered Magnetization B
        
    this->Lnem = 0.0;						// nematic order SA * SB
    this->Ldel = 0.0;						// DA - DB

	this->SA = std::vector<int>(N);
    this->SB = std::vector<int>(N);
    this->DA = std::vector<float>(N);
    this->DB = std::vector<float>(N);
};
	
void SystemState::print_E(){
		std::cout << this->E << std::endl;
};

}