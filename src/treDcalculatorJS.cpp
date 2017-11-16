#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <iostream>
#include <string>
//#include <regex>
#include <algorithm>

std::string exec(const char* cmd);
double samDregion(std::string s_Bamfile, std::string s_PeakCoordinates);

void treDcalculatorJS(double transferArrJS[][1200], std::string s_bamPath, std::string s_geneinfo, double d_minDepth, double depthArr[], int j) {

	std::string s_chrCommand = std::string("printf '") + s_geneinfo + std::string("' | cut -f 1");
	std::string s_startCommand = std::string("printf '") + s_geneinfo + std::string("' | cut -f 2");
	std::string s_endCommand = std::string("printf '") + s_geneinfo + std::string("' | cut -f 3");
	std::string s_strandCommand = std::string("printf '") + s_geneinfo + std::string("' | cut -f 4");

	std::string s_chr = exec(s_chrCommand.c_str());
	std::string s_start = exec(s_startCommand.c_str());
	std::string s_end = exec(s_endCommand.c_str());
	std::string s_strand = exec(s_strandCommand.c_str());

	// remove newline at end of line that find gives
	s_chr.erase(std::remove(s_chr.begin(), s_chr.end(), '\n'), s_chr.end());
	s_start.erase(std::remove(s_start.begin(), s_start.end(), '\n'), s_start.end());
	s_end.erase(std::remove(s_end.begin(), s_end.end(), '\n'), s_end.end());
	s_strand.erase(std::remove(s_strand.begin(), s_strand.end(), '\n'), s_strand.end());


	double d_genelength = atof(s_end.c_str()) - atof(s_start.c_str()); 
	double d_percent = d_genelength/800;

	int i_regionstart;
	int i_regionend;
	int i_upstart;
	int i_upend;
	int i_downstart;
	int i_downend;
	int i_flank = int(d_genelength/4);

	double geneArr[800] = {0}; 	// lowest index is TSS
	double upArr[200] = {0}; 	// lowest index is farthest from gene
	double downArr[200] = {0}; 	// lowest index is TES

	std::string s_geneCoordinates;
	std::string s_upCoordinates;
	std::string s_downCoordinates;
	
	double d_depth = 0;

	for (int i = 0; i < 200; i++) { // upstream
		if (s_strand == "-") {
			i_upstart = atoi(s_end.c_str()) + i_flank - (i_flank/200*i);
			i_upend = atoi(s_end.c_str()) + i_flank - (i_flank/200*(i+1));
			s_upCoordinates = s_chr + ":" + std::to_string(i_upend) + "-" + std::to_string(i_upstart);

			// CALCULATE DEPTH
			d_depth = (samDregion(s_bamPath, s_upCoordinates) / depthArr[j] * d_minDepth);
		}
		if (s_strand == "+") {
			i_upstart = atoi(s_start.c_str()) - i_flank + (i_flank/200*i);
			i_upend = atoi(s_start.c_str()) - i_flank + (i_flank/200*(i+1));
			s_upCoordinates = s_chr + ":" + std::to_string(i_upstart) + "-" + std::to_string(i_upend);

			// CALCULATE DEPTH
			d_depth = (samDregion(s_bamPath, s_upCoordinates) / depthArr[j] * d_minDepth);
		}
		upArr[i] = d_depth;
	}

	for (int i = 0; i < 200; i++) { // downstream
		if (s_strand == "-") {
			i_downstart = atoi(s_start.c_str()) - (i_flank/200*i);
			i_downend = atoi(s_start.c_str()) - (i_flank/200*(i+1));
			s_downCoordinates = s_chr + ":" + std::to_string(i_downend) + "-" + std::to_string(i_downstart);

			// CALCULATE DEPTH
			d_depth = (samDregion(s_bamPath, s_downCoordinates) / depthArr[j] * d_minDepth);
		}
		if (s_strand == "+") {
			i_downstart = atoi(s_end.c_str()) + (i_flank/200*i);
			i_downend = atoi(s_end.c_str()) + (i_flank/200*(i+1));
			s_downCoordinates = s_chr + ":" + std::to_string(i_downstart) + "-" + std::to_string(i_downend);

			// CALCULATE DEPTH
			d_depth = (samDregion(s_bamPath, s_downCoordinates) / depthArr[j] * d_minDepth);
		}
		downArr[i] = d_depth;
	}

	for (int i = 0; i < 800; i++) { // gene body
		if (s_strand == "-") {
			i_regionstart = atoi(s_end.c_str()) - i*d_percent;
			i_regionend = atoi(s_end.c_str()) - (i+1)*d_percent;
			s_geneCoordinates = s_chr + ":" + std::to_string(i_regionend) + "-" + std::to_string(i_regionstart);

			// CALCULATE DEPTH
			d_depth = (samDregion(s_bamPath, s_geneCoordinates) / depthArr[j] * d_minDepth);
		}	
		if (s_strand == "+") {
			i_regionstart = atoi(s_start.c_str()) + i*d_percent;
			i_regionend = atoi(s_start.c_str()) + (i+1)*d_percent;
			s_geneCoordinates = s_chr + ":" + std::to_string(i_regionstart) + "-" + std::to_string(i_regionend);

			// CALCULATE DEPTH
			d_depth = (samDregion(s_bamPath, s_geneCoordinates) / depthArr[j] * d_minDepth);			
		}
		geneArr[i] = d_depth;
	}

	// write to normGeneVector
	for (int i = 0; i < 200; i++) // upstream
		transferArrJS[j][i] = upArr[i]; 

	for (int i = 0; i < 800; i++) // gene body
		transferArrJS[j][i+200] = geneArr[i]; 

	for (int i = 0; i < 200; i++) // downstream
		transferArrJS[j][i+1000] = downArr[i]; 

}
