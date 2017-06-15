#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

std::string exec(const char* cmd);

// Returns a samtools depth double from a specified BAM file 
double samDtotal(std::string s_Bamfile)
{
	std::string s_SamDepth ("samtools depth ");
	std::string s_Awk (" | awk '{sum+=$3} END {print sum}'");		
	std::string s_depth = s_SamDepth + s_Bamfile + s_Awk;

	// Run samtools depth and print to temp file
	// XXX check if s_bamfile has an index file XXX

	std::string s_line;
	try {
		s_line = exec(s_depth.c_str());
	}
	catch(...){
		std::cerr << "Cannot access samtools or bedfile not in standard format" << std::endl;
		std::exit(0);
	}

	double d_line = atof(s_line.c_str());

	return d_line;
}
