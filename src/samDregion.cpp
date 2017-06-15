#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>


std::string exec(const char* cmd);

// Returns a samtools depth double from a specified genome region and BAM file 
double samDregion(std::string s_Bamfile, std::string s_PeakCoordinates)
{
	std::string s_SamDepth ("samtools depth ");
	std::string s_RegionsFlag (" -r ");
	std::string s_Awk (" | awk '{sum+=$3} END {print sum}'");		
	std::string s_PeakDepth = s_SamDepth + s_Bamfile + s_RegionsFlag + s_PeakCoordinates + s_Awk;

	// Run samtools depth at a genome region and print to temp file

	std::string s_line;
	try {
		s_line = exec(s_PeakDepth.c_str());
	}
	catch(...){
		std::cerr << "Cannot access samtools" << std::endl;
		std::exit(0);
	}

	double d_line = atof(s_line.c_str());

	return d_line;
}
