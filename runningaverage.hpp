/**
 * @file runningaverage.hpp
 *
 * defines and prints the variables that are averaged during the simulation
 */
#ifndef RUNNINGAVERAGE_HPP
#define RUNNINGAVERAGE_HPP
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
using namespace std;
namespace runningaverage {
namespace EISH {
	/*
	 Store average energy, magnetization, staggered magn., spin-correlator, ...
	 */
	struct RunningAverage {
		//constructor
		RunningAverage();
		
		string header_output_string();
		string output_string();
		
		void reset();
		
		void increment(double e, double m, double ma, double mb, double sma, double smb, double da, double db, double sda, double sdb, double mamb, double dadb, double mada, double madb, double mbda, double mbdb, double mambda, double mambdb, double madadb, double mbdadb, double mambdadb);
		
		void div(int subpasses);
		
		double avE;
		double avE2;
		double avM;
		double avMA;
		double avMB;
		double avSMA;
		double avSMB;
		double avDA;
		double avDB;
		double avSDA;
		double avSDB;
		double avMAMB;
		double avDADB;
		double avMADA;
		double avMADB;
		double avMBDA;
		double avMBDB;
		double avMAMBDA;
		double avMAMBDB;
		double avMADADB;
		double avMBDADB;
		double avMAMBDADB;
	};
}
}
#endif
