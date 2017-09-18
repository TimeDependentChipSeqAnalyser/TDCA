#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
using namespace std;


void rSlopeRes(std::string s_rPltsName, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, int i_bamFiles, std::string s_name) 
{

	// write r script
	ofstream rLinResResPlt;
	rLinResResPlt.open(s_rPltsName, std::ios::app);

	// write csv file
	ofstream csvLinRes;
	std::string s_LinResCsvName;
	
	s_LinResCsvName = "./" + s_rPltsName + "_data/linRes.csv";

	csvLinRes.open(s_LinResCsvName);

	rLinResResPlt << "linRes <- read.csv(file=\"./" << s_rPltsName << "_data/linRes.csv\", header=FALSE)\n";


	// Obtain residuals from linear regression output file.
	std::string s_txtLinRegName = s_name + ".tdca.LinReg.txt";
	std::string s_line;
	std::ifstream inFile (s_txtLinRegName);
	if (inFile.is_open()) {
		std::getline (inFile, s_line); // header
		for (int j = 0; j < i_peakNumber-1; j++) { 
			std::getline (inFile, s_line);
			std::stringstream linestream(s_line);
			std::string s_linRes;
			for (int j = 0; j < i_bamFiles+5; j++) { 
				std::getline(linestream, s_linRes, '\t');
			}
			std::getline(linestream, s_linRes, '\t');
			csvLinRes << s_linRes << "\n";
		}
	} else { std::cout << "Unable to open file linear regression output file. Program terminated." << std::endl; }

	// last line
	std::getline (inFile, s_line);
	std::stringstream linestream(s_line);
	std::string s_linRes;
	for (int j = 0; j < i_bamFiles+5; j++) { 
		std::getline(linestream, s_linRes, '\t');
	}
	std::getline(linestream, s_linRes, '\t');
	csvLinRes << s_linRes;
	inFile.close();




	
	rLinResResPlt << "dfLinRes<- data.frame(linRes,SlopeBox)\n"; 
	rLinResResPlt << "linRes <- ggplot(dfLinRes, aes(dfLinRes[,1],dfLinRes[,2])) +\n"; 
	rLinResResPlt << "	geom_point(size = 0.3, alpha = 0.25, na.rm = T) +\n";
	rLinResResPlt << "	ggtitle(\"Slope and Residuals from Linear Regression\") +\n";
	rLinResResPlt << "	labs(y = \"Slope from Linear Regression\",size = 14) +\n";
	rLinResResPlt << "	labs(x = \"Residuals from Linear Regression\",size = 14) +\n";


	rLinResResPlt << "	theme(\n";
	rLinResResPlt << "		plot.title = element_text(size=10),\n";
	rLinResResPlt << "		axis.text = element_text(size = 8, colour = \"black\"),\n";
	rLinResResPlt << "		legend.key = element_rect(fill = \"white\"),\n";
	rLinResResPlt << "		legend.background = element_rect(fill = \"white\"),\n";
	rLinResResPlt << "		legend.position = (\"bottom\"),\n";
	rLinResResPlt << "		legend.text = element_text(size = 7),\n";
	rLinResResPlt << "		legend.title=element_blank(),\n";
	rLinResResPlt << "		panel.grid.major = element_line(colour = \"white\"),\n";
	rLinResResPlt << "		panel.grid.minor = element_blank(),\n";
	rLinResResPlt << "		panel.background = element_rect(fill = \"white\"),\n";
	rLinResResPlt << "		axis.line.x = element_line(colour = \"black\"),\n";
	rLinResResPlt << "		axis.line.y = element_line(colour = \"black\")\n";
	rLinResResPlt << "	)\n";


	rLinResResPlt.close();
	csvLinRes.close();
}




