#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

double samDregion(std::string s_Bamfile, std::string s_PeakCoordinates);
std::string exec(const char* cmd);

// Returns a samtools depth double from a specified genome region and BAM file 
double samDnonPeaks(std::string s_Bamfile, std::string s_bed, int i_peakNumber)
{
	std::ifstream bedIn(s_bed);
	std::string s_line;

	// holds loci from variable bed file as such: CHR:START-END
	std::vector<std::string> s_loci;

	while(std::getline(bedIn, s_line)) {  
		std::istringstream iss(s_line);
		std::string s_subLine;
		std::string s_temp;
		std::getline(iss, s_subLine, '\t');
		s_temp = s_subLine + ":";
		std::getline(iss, s_subLine, '\t');
		s_temp += s_subLine + "-";
		std::getline(iss, s_subLine, '\t');
		s_temp += s_subLine;
		s_loci.push_back(s_temp);
	}

	double d_nonPeakDepth = 0;

	for (int i = 0; i < i_peakNumber; i++) {
		d_nonPeakDepth += (samDregion(s_Bamfile, s_loci[i]));
	}
	return d_nonPeakDepth;
}
