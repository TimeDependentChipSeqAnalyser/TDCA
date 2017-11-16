#include <fstream>
#include <iostream>
#include <string>
#include <vector>
//#include <regex>
#include <algorithm>
#include <cmath>
#include <math.h> 
using namespace std;

std::string exec(const char* cmd);

void inflectionWriter(std::vector<std::string> &drcInflectionVector, std::vector<std::string> &drcResidualsVector, std::vector<std::string> &drcUpperAsymVector, std::vector<std::string> &linResVector) 
{
	std::string s_drcScript = "Rscript drcInflection.R > drcInflections.txt";
	std::string s_linearScript = "Rscript linearResiduals.R > linearResiduals.txt";

	try {
		std::system(s_drcScript.c_str());
		std::system(s_linearScript.c_str());
	}
	catch(...){
		std::cerr << "Cannot access R" << std::endl;
		std::exit(0);
	}


	// sigmoidal inflections
	std::remove("inflections.txt");
	std::string s_inflectionsScript = "grep -o 'e:(Intercept)[[:space:]]*[^[:space:]]*' drcInflections.txt | grep -o '[[:space:]]*[^[:space:]]*' | grep -o '[^[:space:]]*' > inflections.txt";

	try {
		std::system(s_inflectionsScript.c_str());
		std::system(s_inflectionsScript.c_str());
	} catch(...){
		std::cerr << "s_inflectionsScript not working." << std::endl;
		std::exit(0);
	}

	std::ifstream inflecFile("inflections.txt");	
	int i_matchCount = 0;
	std::string s_line;	
			
	while (std::getline(inflecFile, s_line)) {
		std::getline(inflecFile, s_line);
		s_line.erase(std::remove(s_line.begin(), s_line.end(), '\n'), s_line.end());
		try {
			std::cout << std::fixed;
			double d = std::stod(s_line);
			if ( (std::isinf(d)) || (std::isnan(d)) ) {
				drcInflectionVector[i_matchCount] = std::string("NaN");
			} else {
				drcInflectionVector[i_matchCount] = std::to_string(d);
			}
		} catch(...){
			drcInflectionVector[i_matchCount] = std::string("NaN");
		}
		i_matchCount++;
	}
	std::remove("inflections.txt");



	// get upper asymptote
	i_matchCount = 0;
	std::remove("asymptotes.txt");
	std::string s_asymptoteScript = "grep -o 'd:(Intercept)[[:space:]]*[^[:space:]]*' drcInflections.txt | grep -o '[[:space:]]*[^[:space:]]*' | grep -o '[^[:space:]]*' > asymptotes.txt";

	try {
		std::system(s_asymptoteScript.c_str());
		std::system(s_asymptoteScript.c_str());
	} catch(...){
		std::cerr << "s_asymptoteScript not working." << std::endl;
		std::exit(0);
	}

	std::ifstream assFile("asymptotes.txt");	

	while (std::getline(assFile, s_line)) {
		std::getline(assFile, s_line);
		s_line.erase(std::remove(s_line.begin(), s_line.end(), '\n'), s_line.end());
		try {
			std::cout << std::fixed;
			double d = std::stod(s_line);
			if ( (std::isinf(d)) || (std::isnan(d)) ) {
				drcUpperAsymVector[i_matchCount] = std::string("NaN");
			} else {
				drcUpperAsymVector[i_matchCount] = std::to_string(d);
			}
		} catch(...){
			drcUpperAsymVector[i_matchCount] = std::string("NaN");
		}
		i_matchCount++;

	}
	std::remove("asymptotes.txt");





	// sygmoidal residuals
	i_matchCount = 0;
	std::remove("asymres.txt");
	std::string s_arScript = "grep 'degrees of freedom' drcInflections.txt | grep -o '[^[:space:]]*[[:space:]]*(' | grep -o '[^[:space:]]*' > asymres.txt";

	try {
		std::system(s_arScript.c_str());
		std::system(s_arScript.c_str());
	} catch(...){
		std::cerr << "s_arScript not working." << std::endl;
		std::exit(0);
	}

	std::ifstream arFile("asymres.txt");	

	while (std::getline(arFile, s_line)) {
		s_line.erase(std::remove(s_line.begin(), s_line.end(), '\n'), s_line.end());
		try {
			std::cout << std::fixed;
			double d = std::stod(s_line);
			if ( (std::isinf(d)) || (std::isnan(d)) ) {
				drcResidualsVector[i_matchCount] = std::string("NaN");
			} else {
				drcResidualsVector[i_matchCount] = std::to_string(d);
			}
		} catch(...){
			drcResidualsVector[i_matchCount] = std::string("NaN");
		}
		std::getline(arFile, s_line);
		i_matchCount++;
	}
	std::remove("asymres.txt");





	// linear residuals
	i_matchCount = 0;
	std::remove("linres.txt");
	std::string s_lrScript = "grep 'Residual standard error:' linearResiduals.txt | grep -o '[^[:space:]]*[[:space:]]*on' | grep -o '[^[:space:]]*' > linres.txt";

	try {
		std::system(s_lrScript.c_str());
		std::system(s_lrScript.c_str());
	} catch(...){
		std::cerr << "s_lrScript not working." << std::endl;
		std::exit(0);
	}
	std::ifstream lrFile("linres.txt");	
	while (std::getline(lrFile, s_line)) {
		s_line.erase(std::remove(s_line.begin(), s_line.end(), '\n'), s_line.end());

		try {
			std::cout << std::fixed;
			double d = std::stod(s_line);
			if ( (std::isinf(d)) || (std::isnan(d)) ) {
				linResVector[i_matchCount] = std::string("NaN");
			} else {
				linResVector[i_matchCount] = std::to_string(d);
			}
		} catch(...){
			linResVector[i_matchCount] = std::string("NaN");
		}
		i_matchCount++;
		std::getline(lrFile, s_line);
	}
	std::remove("linres.txt");




	std::remove("drcInflections.txt");
	std::remove("drcInflection.R");

	std::remove("linearResiduals.txt");
	std::remove("linearResiduals.R");



}

