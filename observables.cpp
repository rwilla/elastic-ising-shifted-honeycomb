#include "observables.hpp"
using namespace std;
namespace observables{
namespace EISH {
	SystemState::SystemState(int N){
	    //initialize structure elements
	    
	    this->sysE = 0.0;
		this->sysM = 0.0;
		this->sysMA = 0.0;
		this->sysMB = 0.0;
		this->sysSMA = 0.0;
		this->sysSMB = 0.0;
		this->sysDA = 0.0;
		this->sysDB = 0.0;
		this->sysSDA = 0.0;
		this->sysSDB = 0.0;
		this->sysMAMB = 0.0;
		this->sysDADB = 0.0;
		this->sysMADA = 0.0;
		this->sysMADB = 0.0;
		this->sysMBDA = 0.0;
		this->sysMBDB = 0.0;
		this->sysMAMBDA = 0.0;
		this->sysMAMBDB = 0.0;
		this->sysMADADB = 0.0;
		this->sysMBDADB = 0.0;
		this->sysMAMBDADB = 0.0;
	    
	    this->pass = 0;							// current pass
		
		this->SA = vector<int>(N);
	    this->SB = vector<int>(N);
	    this->DA = vector<float>(N);
	    this->DB = vector<float>(N);
	}
	
	
	
    void SystemState::reset_sys_vars(){
    	this->sysE = 0.0;
		this->sysM = 0.0;
		this->sysMA = 0.0;
		this->sysMB = 0.0;
		this->sysSMA = 0.0;
		this->sysSMB = 0.0;
		this->sysDA = 0.0;
		this->sysDB = 0.0;
		this->sysSDA = 0.0;
		this->sysSDB = 0.0;
		this->sysMAMB = 0.0;
		this->sysDADB = 0.0;
		this->sysMADA = 0.0;
		this->sysMADB = 0.0;
		this->sysMBDA = 0.0;
		this->sysMBDB = 0.0;
		this->sysMAMBDA = 0.0;
		this->sysMAMBDB = 0.0;
		this->sysMADADB = 0.0;
		this->sysMBDADB = 0.0;
		this->sysMAMBDADB = 0.0;
    }
    	
	void SystemState::div_sys_vars(int N){
		this->sysE /= N;
		this->sysM /= N;
		this->sysMA /= N;
		this->sysMB /= N;
		this->sysSMA /= N;
		this->sysSMB /= N;
		this->sysDA /= N;
		this->sysDB /= N;
		this->sysSDA /= N;
		this->sysSDB /= N;
		this->sysMAMB /= N;
		this->sysDADB /= N;
		this->sysMADA /= N;
		this->sysMADB /= N;
		this->sysMBDA /= N;
		this->sysMBDB /= N;
		this->sysMAMBDA /= N;
		this->sysMAMBDB /= N;
		this->sysMADADB /= N;
		this->sysMBDADB /= N;
		this->sysMAMBDADB /= N;
	}
		
}
}