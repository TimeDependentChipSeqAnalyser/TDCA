#include <iostream>
#include <fstream>
#include <iostream>  

std::string exec(const char* cmd);

void pkgCheck() {

	std::ofstream outFile;
	outFile.open ("pkgCheck.R");

	// check if drc is installed
	outFile << "is.installed <- function(mypkg){\n";
	outFile << "	is.element(mypkg, installed.packages()[,1])\n";
	outFile << "} \n";
	outFile << "if (is.installed(\"drc\")){\n";
	outFile << "	cat(\"TRUE\")\n";
	outFile << "}\n";
	outFile.close();

	std::string s_installed;

	try {
		s_installed = exec(std::string("Rscript pkgCheck.R").c_str());
	}
	catch(...){
		std::cerr << "Cannot access R. Program terminated" << std::endl;
		std::exit(0);
	}

	if (s_installed != "TRUE") {			
		std::cerr << "The R package drc is not installed. Program terminated." << std::endl;
		exec(std::string("rm pkgCheck.R").c_str());
		std::exit(0); 
	}

	exec(std::string("rm pkgCheck.R").c_str());
}
