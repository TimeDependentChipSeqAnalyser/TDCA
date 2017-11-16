#include <fstream>
#include <iostream>
#include <string>
#include <vector>
//#include <regex>
#include <algorithm>
#include <cmath>
#include <math.h> 
#include <sstream>
#include <omp.h>
using namespace std;

std::string exec(const char* cmd);
/**
 argument vector variable recap. This function uses regex to get relevent variables.

 s_drcFI_vec	forward inflection vector
 s_drcFH_vec	forward hill's slope vector
 s_drcFU_vec	forward upper asymptote vector
 s_drcFL_vec	forward lower asymptote vector
 s_drcRI_vec	reverse inflection vector
 s_drcRH_vec	reverse hill's slope vector
 s_drcRU_vec	reverse upper asymptote vector
 s_drcRL_vec	reverse lower asymptote vector
 s_type_vec 	vector holding model type
**/

void inflectionWriterVersitile(int i_maxThreads, int i_peakNumber, std::vector<std::string> &s_drcFI_vec, std::vector<std::string> &s_drcFH_vec, std::vector<std::string> &s_drcFU_vec, std::vector<std::string> &s_drcFL_vec, std::vector<std::string> &s_drcFE_vec, std::vector<std::string> &s_drcRI_vec, std::vector<std::string> &s_drcRH_vec, std::vector<std::string> &s_drcRU_vec, std::vector<std::string> &s_drcRL_vec, std::vector<std::string> &s_drcRE_vec, std::vector<std::string> &s_type_vec, std::vector<std::string> &s_drcRiseRes_vec, std::vector<std::string> &s_drcFallRes_vec, std::string s_name, std::vector<std::vector<std::string> > &transferArray)
{
	if ( i_maxThreads == 1) { std::cout << "Running drc R script." << std::endl; }
	else { std::cout << "Running drc R script using " << i_maxThreads << " processors." << std::endl; }

	try {
		if ( i_maxThreads == 1) {
			std::string s_drcScript = "Rscript drcVersitile." + s_name + ".R > drcVersitile." + s_name + ".txt";
			std::system(s_drcScript.c_str());
		} else {
			#pragma omp parallel for
			for (int i = 1; i < i_maxThreads+1; i++) {
				std::string s_drcScript = "Rscript drcVersitile." + s_name + to_string(i) + ".R > drcVersitile." + s_name + to_string(i) + ".txt";
				std::system(s_drcScript.c_str());
			}
		}
	} catch(...){
		std::cerr << "Cannot access R" << std::endl;
		std::exit(0);
	}

	if ( i_maxThreads > 1) {
		std::string s_cat = "cat ";
		for (int i = 1; i < i_maxThreads+1; i++) { s_cat += "drcVersitile." + s_name + to_string(i) + ".txt "; }
		s_cat += "> drcVersitile." + s_name + ".txt";
		std::system(s_cat.c_str());
		s_cat = "cat ";
		for (int i = 1; i < i_maxThreads+1; i++) { s_cat += "drcVersitile." + s_name + to_string(i) + ".R "; }
		s_cat += "> drcVersitile." + s_name + ".R";
		std::system(s_cat.c_str());
		for (int i = 1; i < i_maxThreads+1; i++) { 
			std::string s_rm = "rm drcVersitile." + s_name + to_string(i) + ".txt; rm drcVersitile." + s_name + to_string(i) + ".R"; 
			std::system(s_rm.c_str());
		}
	}



	// write parameters to a vector
	std::remove(std::string("parameter." + s_name + ".txt").c_str());
	std::string s_parameter;

	// for progress
	float f_progress = 0.0;
	int i_barWidth = 20;
	double d_progiterator = 0;

	for (int i = 0; i < 5; i++) { 
		
		if (i == 0) { s_parameter = 'b';} // hill's slope
		if (i == 1) { s_parameter = 'c';} // lower asymptote
		if (i == 2) { s_parameter = 'd';} // upper asymptote
		if (i == 3) { s_parameter = 'e';} // inflection point
		if (i == 4) { s_parameter = 'f';} // symetricity value

		std::string s_parameterScript = "grep -o '" + s_parameter + ":(Intercept)[[:space:]]*[^[:space:]]*'" +
			" drcVersitile." + s_name + ".txt | grep -o '[[:space:]]*[^[:space:]]*' | grep -o '[^[:space:]]*' > parameter." + s_name + ".txt";

		try {
			std::system(s_parameterScript.c_str());
			std::system(s_parameterScript.c_str());
		} catch(...){
			std::cerr << "s_parameterScript failed on ";
			if (i == 0) { std::cerr << "hill's coefficient calculation." << std::endl; }
			if (i == 1) { std::cerr << "lower asymptote calculation." << std::endl; }
			if (i == 2) { std::cerr << "upper asymptote calculation." << std::endl; }
			if (i == 3) { std::cerr << "inflection point calculation." << std::endl; }
			if (i == 4) { std::cerr << "symetricity calculation." << std::endl; }
			std::exit(0);
		}

		std::ifstream inFile("parameter." + s_name + ".txt");	
		int i_matchCount = 0;

		std::string s_line;	
		for (int j = i_matchCount; j < i_peakNumber; j++) {
			std::getline(inFile, s_line);
			std::getline(inFile, s_line);
			s_line.erase(std::remove(s_line.begin(), s_line.end(), '\n'), s_line.end());
			if (s_type_vec[i_matchCount] == "rise") { // XXX HERE
				try {
					std::cout << std::fixed;
					double d = std::stod(s_line);
					if ( (std::isinf(d)) || (std::isnan(d)) ) {
						// ALTER VECTOR ACCORDINGLY
						if (i == 0) { s_drcFH_vec[i_matchCount] = std::string("NaN"); }
						if (i == 1) { s_drcFL_vec[i_matchCount] = std::string("NaN"); }
						if (i == 2) { s_drcFU_vec[i_matchCount] = std::string("NaN"); }
						if (i == 3) { s_drcFI_vec[i_matchCount] = std::string("NaN"); }
						if (i == 4) { s_drcFE_vec[i_matchCount] = std::string("NaN"); }
					} else {
						if (i == 0) { s_drcFH_vec[i_matchCount] = std::to_string(d); }
						if (i == 1) { s_drcFL_vec[i_matchCount] = std::to_string(d); }
						if (i == 2) { s_drcFU_vec[i_matchCount] = std::to_string(d); }
						if (i == 3) { s_drcFI_vec[i_matchCount] = std::to_string(d); }
						if (i == 4) { s_drcFE_vec[i_matchCount] = std::to_string(d); }
					} 
				} catch(...) {
					if (i == 0) { s_drcFH_vec[i_matchCount] = std::string("NaN"); }
					if (i == 1) { s_drcFL_vec[i_matchCount] = std::string("NaN"); }
					if (i == 2) { s_drcFU_vec[i_matchCount] = std::string("NaN"); }
					if (i == 3) { s_drcFI_vec[i_matchCount] = std::string("NaN"); }
					if (i == 4) { s_drcFE_vec[i_matchCount] = std::string("NaN"); }
				} //try
				if ( (i == 1) && (s_line == "FORCED\"") ) { s_drcFL_vec[i_matchCount] = std::string("forced to 0"); }
			} else if (s_type_vec[i_matchCount] == "fall") { 
				try {
					std::cout << std::fixed;
					double d = std::stod(s_line);
					if ( (std::isinf(d)) || (std::isnan(d)) ) {
						// ALTER VECTOR ACCORDINGLY
						if (i == 0) { s_drcRH_vec[i_matchCount] = std::string("NaN"); }
						if (i == 1) { s_drcRL_vec[i_matchCount] = std::string("NaN"); }
						if (i == 2) { s_drcRU_vec[i_matchCount] = std::string("NaN"); }
						if (i == 3) { s_drcRI_vec[i_matchCount] = std::string("NaN"); }
						if (i == 4) { s_drcRE_vec[i_matchCount] = std::string("NaN"); }
					} else {
						if (i == 0) { s_drcRH_vec[i_matchCount] = std::to_string(d); }
						if (i == 1) { s_drcRL_vec[i_matchCount] = std::to_string(d); }
						if (i == 2) { s_drcRU_vec[i_matchCount] = std::to_string(d); }
						if (i == 3) { s_drcRI_vec[i_matchCount] = std::to_string(d); }
						if (i == 4) { s_drcRE_vec[i_matchCount] = std::to_string(d); }
					} 
				} catch(...) {
					if (i == 0) { s_drcRH_vec[i_matchCount] = std::string("NaN"); }
					if (i == 1) { s_drcRL_vec[i_matchCount] = std::string("NaN"); }
					if (i == 2) { s_drcRU_vec[i_matchCount] = std::string("NaN"); }
					if (i == 3) { s_drcRI_vec[i_matchCount] = std::string("NaN"); }
					if (i == 4) { s_drcRE_vec[i_matchCount] = std::string("NaN"); }
				} //try
				if ( (i == 1) && (s_line == "FORCED\"") ) { s_drcRL_vec[i_matchCount] = std::string("forced to 0"); }
			} else if (s_type_vec[i_matchCount] == "hill") { 
				try {
					std::cout << std::fixed;
					double d = std::stod(s_line);
					if ( (std::isinf(d)) || (std::isnan(d)) ) {
						// ALTER VECTOR ACCORDINGLY
						if (i == 0) { s_drcFH_vec[i_matchCount] = std::string("NaN"); }
						if (i == 1) { s_drcFL_vec[i_matchCount] = std::string("NaN"); }
						if (i == 2) { s_drcFU_vec[i_matchCount] = std::string("NaN"); }
						if (i == 3) { s_drcFI_vec[i_matchCount] = std::string("NaN"); }
						if (i == 4) { s_drcFE_vec[i_matchCount] = std::string("NaN"); }
					} else {
						if (i == 0) { s_drcFH_vec[i_matchCount] = std::to_string(d); }
						if (i == 1) { s_drcFL_vec[i_matchCount] = std::to_string(d); }
						if (i == 2) { s_drcFU_vec[i_matchCount] = std::to_string(d); }
						if (i == 3) { s_drcFI_vec[i_matchCount] = std::to_string(d); }
						if (i == 4) { s_drcFE_vec[i_matchCount] = std::to_string(d); }
					} 
				} catch(...) {
					if (i == 0) { s_drcFH_vec[i_matchCount] = std::string("NaN"); }
					if (i == 1) { s_drcFL_vec[i_matchCount] = std::string("NaN"); }
					if (i == 2) { s_drcFU_vec[i_matchCount] = std::string("NaN"); }
					if (i == 3) { s_drcFI_vec[i_matchCount] = std::string("NaN"); }
					if (i == 4) { s_drcFE_vec[i_matchCount] = std::string("NaN"); }
				} //try
				if ( (i == 1) && (s_line == "FORCED\"") ) { s_drcFL_vec[i_matchCount] = std::string("forced to 0"); }
				std::getline(inFile, s_line);
				std::getline(inFile, s_line);
				s_line.erase(std::remove(s_line.begin(), s_line.end(), '\n'), s_line.end());
				try {
					std::cout << std::fixed;
					double d = std::stod(s_line);
					if ( (std::isinf(d)) || (std::isnan(d)) ) {
						// ALTER VECTOR ACCORDINGLY
						if (i == 0) { s_drcRH_vec[i_matchCount] = std::string("NaN"); }
						if (i == 1) { s_drcRL_vec[i_matchCount] = std::string("NaN"); }
						if (i == 2) { s_drcRU_vec[i_matchCount] = std::string("NaN"); }
						if (i == 3) { s_drcRI_vec[i_matchCount] = std::string("NaN"); }
						if (i == 4) { s_drcRE_vec[i_matchCount] = std::string("NaN"); }
					} else {
						if (i == 0) { s_drcRH_vec[i_matchCount] = std::to_string(d); }
						if (i == 1) { s_drcRL_vec[i_matchCount] = std::to_string(d); }
						if (i == 2) { s_drcRU_vec[i_matchCount] = std::to_string(d); }
						if (i == 3) { s_drcRI_vec[i_matchCount] = std::to_string(d); }
						if (i == 4) { s_drcRE_vec[i_matchCount] = std::to_string(d); }
					} 
				} catch(...) {
					if (i == 0) { s_drcRH_vec[i_matchCount] = std::string("NaN"); }
					if (i == 1) { s_drcRL_vec[i_matchCount] = std::string("NaN"); }
					if (i == 2) { s_drcRU_vec[i_matchCount] = std::string("NaN"); }
					if (i == 3) { s_drcRI_vec[i_matchCount] = std::string("NaN"); }
					if (i == 4) { s_drcRE_vec[i_matchCount] = std::string("NaN"); }
				} //try
				if ( (i == 1) && (s_line == "FORCED\"") ) { s_drcRL_vec[i_matchCount] = std::string("forced to 0"); }
			} else if (s_type_vec[i_matchCount] == "valley") { 
				try {
					std::cout << std::fixed;
					double d = std::stod(s_line);
					if ( (std::isinf(d)) || (std::isnan(d)) ) {
						// ALTER VECTOR ACCORDINGLY
						if (i == 0) { s_drcRH_vec[i_matchCount] = std::string("NaN"); }
						if (i == 1) { s_drcRL_vec[i_matchCount] = std::string("NaN"); }
						if (i == 2) { s_drcRU_vec[i_matchCount] = std::string("NaN"); }
						if (i == 3) { s_drcRI_vec[i_matchCount] = std::string("NaN"); }
						if (i == 4) { s_drcRE_vec[i_matchCount] = std::string("NaN"); }
					} else {
						if (i == 0) { s_drcRH_vec[i_matchCount] = std::to_string(d); }
						if (i == 1) { s_drcRL_vec[i_matchCount] = std::to_string(d); }
						if (i == 2) { s_drcRU_vec[i_matchCount] = std::to_string(d); }
						if (i == 3) { s_drcRI_vec[i_matchCount] = std::to_string(d); }
						if (i == 4) { s_drcRE_vec[i_matchCount] = std::to_string(d); }
					} 
				} catch(...) {
					if (i == 0) { s_drcRH_vec[i_matchCount] = std::string("NaN"); }
					if (i == 1) { s_drcRL_vec[i_matchCount] = std::string("NaN"); }
					if (i == 2) { s_drcRU_vec[i_matchCount] = std::string("NaN"); }
					if (i == 3) { s_drcRI_vec[i_matchCount] = std::string("NaN"); }
					if (i == 4) { s_drcRE_vec[i_matchCount] = std::string("NaN"); }
				} //try
				if ( (i == 1) && (s_line == "FORCED\"") ) { s_drcRL_vec[i_matchCount] = std::string("forced to 0"); }
				std::getline(inFile, s_line);
				std::getline(inFile, s_line);
				s_line.erase(std::remove(s_line.begin(), s_line.end(), '\n'), s_line.end());
				try {
					std::cout << std::fixed;
					double d = std::stod(s_line);
					if ( (std::isinf(d)) || (std::isnan(d)) ) {
						// ALTER VECTOR ACCORDINGLY
						if (i == 0) { s_drcFH_vec[i_matchCount] = std::string("NaN"); }
						if (i == 1) { s_drcFL_vec[i_matchCount] = std::string("NaN"); }
						if (i == 2) { s_drcFU_vec[i_matchCount] = std::string("NaN"); }
						if (i == 3) { s_drcFI_vec[i_matchCount] = std::string("NaN"); }
						if (i == 4) { s_drcFE_vec[i_matchCount] = std::string("NaN"); }
					} else {
						if (i == 0) { s_drcFH_vec[i_matchCount] = std::to_string(d); }
						if (i == 1) { s_drcFL_vec[i_matchCount] = std::to_string(d); }
						if (i == 2) { s_drcFU_vec[i_matchCount] = std::to_string(d); }
						if (i == 3) { s_drcFI_vec[i_matchCount] = std::to_string(d); }
						if (i == 4) { s_drcFE_vec[i_matchCount] = std::to_string(d); }
					} 
				} catch(...) {
					if (i == 0) { s_drcFH_vec[i_matchCount] = std::string("NaN"); }
					if (i == 1) { s_drcFL_vec[i_matchCount] = std::string("NaN"); }
					if (i == 2) { s_drcFU_vec[i_matchCount] = std::string("NaN"); }
					if (i == 3) { s_drcFI_vec[i_matchCount] = std::string("NaN"); }
					if (i == 4) { s_drcFE_vec[i_matchCount] = std::string("NaN"); }
				} //try
				if ( (i == 1) && (s_line == "FORCED\"") ) { s_drcFL_vec[i_matchCount] = std::string("forced to 0"); }
			} else { // type is undefined
				try {
					std::cout << std::fixed;
					double d = std::stod(s_line);
					if ( (std::isinf(d)) || (std::isnan(d)) ) {
						// ALTER VECTOR ACCORDINGLY
						if (i == 0) { s_drcFH_vec[i_matchCount] = std::string("NaN"); }
						if (i == 1) { s_drcFL_vec[i_matchCount] = std::string("NaN"); }
						if (i == 2) { s_drcFU_vec[i_matchCount] = std::string("NaN"); }
						if (i == 3) { s_drcFI_vec[i_matchCount] = std::string("NaN"); }
						if (i == 4) { s_drcFE_vec[i_matchCount] = std::string("NaN"); }
					} else { // THESE STATEMENTS DETERMINE WHICH CONTAINER TO WRITE TO (FALL OR RISE)
						if ( (i == 0) && (d < 0) ) { s_drcFH_vec[i_matchCount] = std::to_string(d); } // rise
						if ( (i == 0) && (d > 0) ) { s_drcRH_vec[i_matchCount] = std::to_string(d); } // fall
						if ( (i == 1) && (s_drcFH_vec[i_matchCount] != "-") ) // rise
							s_drcFL_vec[i_matchCount] = std::to_string(d);
						if ( (i == 1) && (s_drcRH_vec[i_matchCount] != "-") ) // fall
							s_drcRL_vec[i_matchCount] = std::to_string(d);
						if ( (i == 2) && (s_drcFH_vec[i_matchCount] != "-") ) // rise
							s_drcFU_vec[i_matchCount] = std::to_string(d);
						if ( (i == 2) && (s_drcRH_vec[i_matchCount] != "-") ) // fall
							s_drcRU_vec[i_matchCount] = std::to_string(d);
						if ( (i == 3) && (s_drcFH_vec[i_matchCount] != "-") ) // rise
							s_drcFI_vec[i_matchCount] = std::to_string(d);
						if ( (i == 3) && (s_drcRH_vec[i_matchCount] != "-") ) // fall
							s_drcRI_vec[i_matchCount] = std::to_string(d);
						if ( (i == 4) && (s_drcFH_vec[i_matchCount] != "-") ) // rise
							s_drcFE_vec[i_matchCount] = std::to_string(d);
						if ( (i == 4) && (s_drcRH_vec[i_matchCount] != "-") ) // fall
							s_drcRE_vec[i_matchCount] = std::to_string(d);
					} 
				} catch(...) {
					if (i == 0) { s_drcFH_vec[i_matchCount] = std::string("NaN"); }
					if ( (i == 1) && (s_drcFH_vec[i_matchCount] != "-") ) { s_drcFL_vec[i_matchCount] = "forced to 0"; }
					if ( (i == 1) && (s_drcRH_vec[i_matchCount] != "-") ) { s_drcRL_vec[i_matchCount] = "forced to 0"; }
					if ( (i == 1) && (s_drcFH_vec[i_matchCount] == "NaN") ) { s_drcFL_vec[i_matchCount] = "NaN"; }
					if ( (i == 1) && (s_drcRH_vec[i_matchCount] == "NaN") ) { s_drcRL_vec[i_matchCount] = "NaN"; }
					if (i == 2) { s_drcFU_vec[i_matchCount] = std::string("NaN"); }
					if (i == 3) { s_drcFI_vec[i_matchCount] = std::string("NaN"); }
					if (i == 4) { s_drcFE_vec[i_matchCount] = std::string("NaN"); }
				} //try
			}
			i_matchCount++;

			// Progress
			std::cout << "[";
			int i_pos = i_barWidth * f_progress;
			for (int k = 0; k < i_barWidth; k++) {
				if (k <= i_pos) std::cout << "=";
				else std::cout << " ";
			}
			if ( (i*j) == ((i_peakNumber-1)*4) ) { // 100% complete
				std::cout << "] " << int(100) << " %\r";
				std::cout.flush();
			} else {
				std::cout << "] " << int((f_progress+0.01) * 100.0) << " %\r";
				std::cout.flush();
			}
			f_progress = d_progiterator/((i_peakNumber)*4);
			d_progiterator++;
		} //for
		std::remove(std::string("parameter." + s_name + ".txt").c_str());

	} //for


	// get the residuals and write to vector:
	std::remove(std::string("residuals." + s_name + ".txt").c_str());
	std::string s_residualsScript = "grep 'degrees of freedom)' drcVersitile." + s_name + ".txt | sed 's/(.*//' | grep -o '[^[:space:]]*' | grep -v ']' | grep -v '\"' > residuals." + s_name + ".txt";
	try {
		std::system(s_residualsScript.c_str());
	} catch(...){
		std::cerr << "s_residualsScript failed.";
		std::exit(0);
	}
	std::ifstream resFile("residuals." + s_name + ".txt");	
	std::string s_line;


	int i_iterator = 0;

	for (int i = 0; i < i_peakNumber; i++) {
		std::getline(resFile, s_line);
		if (s_type_vec[i] == "rise") { // Rise
			i_iterator++;
			s_drcRiseRes_vec[i] = s_line;
		} else if (s_type_vec[i] == "fall") { // Fall
			i_iterator++;
			s_drcFallRes_vec[i] = s_line;
		} else if (s_type_vec[i] == "hill") { // Hill
			i_iterator++;
			i_iterator++;
			s_drcRiseRes_vec[i] = s_line;
			std::getline(resFile, s_line);
			s_drcFallRes_vec[i] = s_line;
		} else if (s_type_vec[i] == "valley") { // Valley
			i_iterator++;
			i_iterator++;
			s_drcFallRes_vec[i] = s_line;
			std::getline(resFile, s_line);
			s_drcRiseRes_vec[i] = s_line;
		} else { // Undef
			i_iterator++;
			if (s_drcFH_vec[i] != "-") { // rise
				s_drcRiseRes_vec[i] = s_line;
			} else { // fall
				s_drcFallRes_vec[i] = s_line;
			}
		}
	}

	std::remove(std::string("residuals." + s_name + ".txt").c_str());


	// initialize vectors that hold standard errors of parameters 
	std::vector<std::string> s_bErrors(i_iterator, "-");
	std::vector<std::string> s_cErrors(i_iterator, "-");
	std::vector<std::string> s_dErrors(i_iterator, "-");
	std::vector<std::string> s_eErrors(i_iterator, "-");
	std::vector<std::string> s_fErrors(i_iterator, "-");

	// get the error for each modelling parameter
	std::remove(std::string("parameterErrors." + s_name + ".txt").c_str());
	for (int i = 0; i < 5; i++) {
		if (i == 0) { s_parameter = 'b';} // hill's slope
		if (i == 1) { s_parameter = 'c';} // lower asymptote
		if (i == 2) { s_parameter = 'd';} // upper asymptote
		if (i == 3) { s_parameter = 'e';} // inflection point
		if (i == 4) { s_parameter = 'f';} // symetricity value

		std::string s_parameterErrorsScript = "grep '" + s_parameter + ":(Intercept)' drcVersitile." + s_name + ".txt | sed -n -e 's/^.*) //p' | sed -e 's/^[ \\t]*//' | sed 's/ \\+ /\\t/g'| tr ' ' \\\\t > parameterErrors." + s_name + ".txt";
		try {
			std::system(s_parameterErrorsScript.c_str());
		} catch(...){
			std::cerr << "s_parameterErrorsScript failed.";
			std::exit(0);
		}

		std::ifstream errorFile("parameterErrors." + s_name + ".txt");
		for (int j = 0; j < i_iterator; j++) {
			std::getline(errorFile, s_line);	
			std::stringstream linestream(s_line);
			std::getline(linestream, s_line, '\t');
			std::getline(linestream, s_line, '\t');
			
			if ( (s_line.size() == 0) || (s_line == "ERROR\"") || (s_line == "FORCED\"") || (s_line == "1\"") ) { // if s_line does not exist (no std. error) or is ERROR"
				if (i == 0) { s_bErrors[j] = "NA";} // hill's slope
				if (i == 1) { s_cErrors[j] = "NA";} // lower asymptote
				if (i == 2) { s_dErrors[j] = "NA";} // upper asymptote
				if (i == 3) { s_eErrors[j] = "NA";} // inflection point
				if (i == 4) { s_fErrors[j] = "NA";} // symetricity value
			} else {
				if (i == 0) { s_bErrors[j] = s_line;} // hill's slope
				if (i == 1) { s_cErrors[j] = s_line;} // lower asymptote
				if (i == 2) { s_dErrors[j] = s_line;} // upper asymptote
				if (i == 3) { s_eErrors[j] = s_line;} // inflection point
				if (i == 4) { s_fErrors[j] = s_line;} // symetricity value
			}
		}
	}

	std::remove(std::string("parameterErrors." + s_name + ".txt").c_str());

	// Print std. error of parameters transferArray
	int i_count = 0;
	for (int i = 0; i < i_peakNumber; i++) {
		if (s_type_vec[i] == "rise") { // Rise
			transferArray[i][0] = s_eErrors[i_count];
			transferArray[i][1] = s_bErrors[i_count];
			transferArray[i][2] = s_dErrors[i_count];
			transferArray[i][3] = s_cErrors[i_count];
			transferArray[i][4] = s_fErrors[i_count];
			transferArray[i][5] = "NA";
			transferArray[i][6] = "NA";
			transferArray[i][7] = "NA";
			transferArray[i][8] = "NA";
			transferArray[i][9] = "NA";
			i_count++;
		} else if (s_type_vec[i] == "fall") { // Fall
			transferArray[i][0] = "NA";
			transferArray[i][1] = "NA";
			transferArray[i][2] = "NA";
			transferArray[i][3] = "NA";
			transferArray[i][4] = "NA";
			transferArray[i][5] = s_eErrors[i_count];
			transferArray[i][6] = s_bErrors[i_count];
			transferArray[i][7] = s_dErrors[i_count];
			transferArray[i][8] = s_cErrors[i_count];
			transferArray[i][9] = s_fErrors[i_count];
			i_count++;
		} else if (s_type_vec[i] == "hill") { // Hill
			transferArray[i][0] = s_eErrors[i_count];
			transferArray[i][1] = s_bErrors[i_count];
			transferArray[i][2] = s_dErrors[i_count];
			transferArray[i][3] = s_cErrors[i_count];
			transferArray[i][4] = s_fErrors[i_count];
			transferArray[i][5] = s_eErrors[i_count+1];
			transferArray[i][6] = s_bErrors[i_count+1];
			transferArray[i][7] = s_dErrors[i_count+1];
			transferArray[i][8] = s_cErrors[i_count+1];
			transferArray[i][9] = s_fErrors[i_count+1];
			i_count++; i_count++;
		} else if (s_type_vec[i] == "valley") { // Valley
			transferArray[i][0] = s_eErrors[i_count+1];
			transferArray[i][1] = s_bErrors[i_count+1];
			transferArray[i][2] = s_dErrors[i_count+1];
			transferArray[i][3] = s_cErrors[i_count+1];
			transferArray[i][4] = s_fErrors[i_count+1];
			transferArray[i][5] = s_eErrors[i_count];
			transferArray[i][6] = s_bErrors[i_count];
			transferArray[i][7] = s_dErrors[i_count];
			transferArray[i][8] = s_cErrors[i_count];
			transferArray[i][9] = s_fErrors[i_count];
			i_count++; i_count++;
		} else { // Undef
			if (s_drcFH_vec[i] != "-") { // rise
				transferArray[i][0] = s_eErrors[i_count];
				transferArray[i][1] = s_bErrors[i_count];
				transferArray[i][2] = s_dErrors[i_count];
				transferArray[i][3] = s_cErrors[i_count];
				transferArray[i][4] = s_fErrors[i_count];
				transferArray[i][5] = "NA";
				transferArray[i][6] = "NA";
				transferArray[i][7] = "NA";
				transferArray[i][8] = "NA";
				transferArray[i][9] = "NA";
			} else { // fall
				transferArray[i][0] = "NA";
				transferArray[i][1] = "NA";
				transferArray[i][2] = "NA";
				transferArray[i][3] = "NA";
				transferArray[i][4] = "NA";
				transferArray[i][5] = s_eErrors[i_count];
				transferArray[i][6] = s_bErrors[i_count];
				transferArray[i][7] = s_dErrors[i_count];
				transferArray[i][8] = s_cErrors[i_count];
				transferArray[i][9] = s_fErrors[i_count];
			}
		i_count++;
		}
	}


	std::cout << std::endl;
	std::remove(std::string("drcVersitile." + s_name + ".txt").c_str()); // COMMENT THIS OUT IF YOU WOULD LIKE TO SEE DRC OUTPUT
	std::remove(std::string("drcVersitile." + s_name + ".R").c_str());   // COMMENT THIS OUT IF YOU WOULD LIKE TO SEE DRC OUTPUT
}
