#include "runningaverage.hpp"
using namespace std;
namespace runningaverage {
namespace EISH {
	RunningAverage::RunningAverage(){
		
		this->avE = 0.0;
		this->avE2 = 0.0;
		this->avM = 0.0;
		this->avMA = 0.0;
		this->avMB = 0.0;
		this->avSMA = 0.0;
		this->avSMB = 0.0;
		this->avDA = 0.0;
		this->avDB = 0.0;
		this->avSDA = 0.0;
		this->avSDB = 0.0;
		this->avMAMB = 0.0;
		this->avDADB = 0.0;
		this->avMADA = 0.0;
		this->avMADB = 0.0;
		this->avMBDA = 0.0;
		this->avMBDB = 0.0;
		this->avMAMBDA = 0.0;
		this->avMAMBDB = 0.0;
		this->avMADADB = 0.0;
		this->avMBDADB = 0.0;
		this->avMAMBDADB = 0.0;
		
	}
	
	string RunningAverage::header_output_string(){
		//string str = "pass \t Temperature \t Hz \t <E> \t <E^2> \t <MA> \t <MA^2> \t <MB> \t <MB^2> \t <SMA> \t <SMA^2> \t <SMB> \t <SMB^2> \t <LocCorr> \t <LocCorr^2> \t <LocDisp> \t <LocDisp^2> \t <LocNem> \t <LocNem^2>";
		string str = "pass \t Temperature \t Hz \t <e> \t <e^2> \t <m> \t <ma> \t <mb> \t <sma> \t <smb> \t <da> \t <db> \t <sda> \t <sdb> \t <mamb> \t <dadb> \t <mada> \t <madb> \t <mbda> \t <mbdb> \t <mambda> \t <mambdb> \t <madadb> \t <mbdadb> \t <mambdadb>";
		return str;
	}

	string RunningAverage::output_string(){
		
		stringstream stream;
		stream << fixed << setprecision(12) << this->avE << "\t" << this->avE2 << "\t" << this->avM << "\t" << \
		       this->avMA << "\t" << this->avMB << "\t" << \
		       this->avSMB << "\t" << this->avSMB << "\t" << \
		       this->avDA << "\t" << this->avDB << "\t" << \
		       this->avSDA << "\t" << this->avSDB << "\t" << \
		       this->avMAMB << "\t" << this->avDADB << "\t" << \
		       this->avMADA << "\t" << this->avMADB << "\t" << \
		       this->avMBDA << "\t" << this->avMBDB << "\t" << \
		       this->avMAMBDA << "\t" << this->avMAMBDB << "\t" << \
		       this->avMADADB << "\t" << this->avMBDADB << "\t" << \
		       this->avMAMBDADB;
		string str = stream.str();       
		
		return str;
	}
		
		
	void RunningAverage::reset(){
		
		this->avE = 0.0;
		this->avE2 = 0.0;
		this->avM = 0.0;
		this->avMA = 0.0;
		this->avMB = 0.0;
		this->avSMA = 0.0;
		this->avSMB = 0.0;
		this->avDA = 0.0;
		this->avDB = 0.0;
		this->avSDA = 0.0;
		this->avSDB = 0.0;
		this->avMAMB = 0.0;
		this->avDADB = 0.0;
		this->avMADA = 0.0;
		this->avMADB = 0.0;
		this->avMBDA = 0.0;
		this->avMBDB = 0.0;
		this->avMAMBDA = 0.0;
		this->avMAMBDB = 0.0;
		this->avMADADB = 0.0;
		this->avMBDADB = 0.0;
		this->avMAMBDADB = 0.0;
	}
	
	void RunningAverage::increment(double e, double m, double ma, double mb, double sma, double smb, double da, double db, double sda, double sdb, double mamb, double dadb, double mada, double madb, double mbda, double mbdb, double mambda, double mambdb, double madadb, double mbdadb, double mambdadb){
		
		this->avE += e;
		this->avE2 += e * e;
		this->avM += m;
		this->avMA += ma;
		this->avMB += mb;
		this->avSMA += sma;
		this->avSMB += smb;
		this->avDA += da;
		this->avDB += db;
		this->avSDA += sda;
		this->avSDB += sdb;
		this->avMAMB += mamb;
		this->avDADB += dadb;
		this->avMADA += mada;
		this->avMADB += madb;
		this->avMBDA += mbda;
		this->avMBDB += mbdb;
		this->avMAMBDA += mambda;
		this->avMAMBDB += mambdb;
		this->avMADADB += madadb;
		this->avMBDADB += mbdadb;
		this->avMAMBDADB += mambdadb;
	}
	
	void RunningAverage::div(int subpasses){
		
		this->avE /= subpasses;
		this->avE2 /= subpasses;
		this->avM /= subpasses;
		this->avMA /= subpasses;
		this->avMB /= subpasses;
		this->avSMA /= subpasses;
		this->avSMB /= subpasses;
		this->avDA /= subpasses;
		this->avDB /= subpasses;
		this->avSDA /= subpasses;
		this->avSDB /= subpasses;
		this->avMAMB /= subpasses;
		this->avDADB /= subpasses;
		this->avMADA /= subpasses;
		this->avMADB /= subpasses;
		this->avMBDA /= subpasses;
		this->avMBDB /= subpasses;
		this->avMAMBDA /= subpasses;
		this->avMAMBDB /= subpasses;
		this->avMADADB /= subpasses;
		this->avMBDADB /= subpasses;
		this->avMAMBDADB /= subpasses;
	}
	
}
}
