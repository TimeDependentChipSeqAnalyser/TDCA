#include <string>
#include <iostream>
#include <vector>
//#include <regex>
#include <fstream> 
#include <algorithm>

std::string exec(const char* cmd);

std::string geneToCoordinates(std::string s_gene, std::string s_genelistFile, std::string s_name)
{

	std::string s_grep = "grep '" + s_gene + "' " + s_genelistFile + " > grep.gene." + s_name + ".txt";

	try {
		std::system(s_grep.c_str());
	}
	catch(...){
		std::cerr << "Grep not functioning." << std::endl;
		std::exit(0);
	}

	std::ifstream inList("grep.gene." + s_name + ".txt"); 
	std::string s_line;

	std::string s_coordinates;

	while (std::getline(inList, s_line, '\n')) {
		
		std::string s_command = std::string("printf '") + s_line + std::string("' | cut -f 2,3,4,5");
		s_coordinates = exec(s_command.c_str());
		s_coordinates.erase(std::remove(s_coordinates.begin(), s_coordinates.end(), '\n'), s_coordinates.end());
		if (s_coordinates != "") {
			std::remove(std::string("grep.gene." + s_name + ".txt").c_str());
			return s_coordinates;	
		}
	}

	std::remove(std::string("grep.gene." + s_name + ".txt").c_str());
}
























