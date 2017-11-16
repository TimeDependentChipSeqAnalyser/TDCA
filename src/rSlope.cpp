#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
using namespace std;


void rSlope(std::string s_rPltsName, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, int i_bamFiles, std::string s_name) 
{

	// write r script
	ofstream rSlopeBoxPlt;
	rSlopeBoxPlt.open(s_rPltsName, std::ios::app);

	// write csv file
	ofstream csvSlope;
	std::string s_SlopeCsvName;
	
	s_SlopeCsvName = "./" + s_rPltsName + "_data/SlopeBox.csv";

	csvSlope.open(s_SlopeCsvName);

	rSlopeBoxPlt << "SlopeBox <- read.csv(file=\"./" << s_rPltsName << "_data/SlopeBox.csv\", header=FALSE)\n";


	// Obtain slope from linear regression output file.
	std::string s_txtLinRegName = s_name + ".tdca.LinReg.txt";
	std::string s_line;
	std::ifstream inFile (s_txtLinRegName);
	if (inFile.is_open()) {
		std::getline (inFile, s_line); // header
		for (int j = 0; j < i_peakNumber-1; j++) { 
			std::getline (inFile, s_line);
			std::stringstream linestream(s_line);
			std::string s_slope;
			for (int j = 0; j < i_bamFiles+4; j++) { 
				std::getline(linestream, s_slope, '\t');
			}
			std::getline(linestream, s_slope, '\t');
			csvSlope << s_slope << "\n";
		}
	} else { std::cout << "Unable to open file linear regression output file. Program terminated." << std::endl; }

	// last line
	std::getline (inFile, s_line);
	std::stringstream linestream(s_line);
	std::string s_slope;
	for (int j = 0; j < i_bamFiles+4; j++) { 
		std::getline(linestream, s_slope, '\t');
	}
	std::getline(linestream, s_slope, '\t');
	csvSlope << s_slope;
	inFile.close();
	
	rSlopeBoxPlt << "dfSlopeBox<- data.frame(idDCov,SlopeBox)\n"; 
	rSlopeBoxPlt << "slope <- ggplot(dfSlopeBox, aes(dfSlopeBox[,1],dfSlopeBox[,2], fill=dfSlopeBox[,1])) +\n"; 
	rSlopeBoxPlt << "	geom_boxplot(outlier.shape = NA) +\n";
	rSlopeBoxPlt << "	scale_fill_brewer(palette=\"Dark2\") +\n";
	rSlopeBoxPlt << "	ggtitle(\"Slope from Linear Regression by Locus Type\") +\n";
	rSlopeBoxPlt << "	labs(y = \"Slope from Linear Regression\",size = 14) +\n";
	rSlopeBoxPlt << "	labs(x = \"Locus Type\",size = 14) +\n";
	rSlopeBoxPlt << "	theme(\n";
	rSlopeBoxPlt << "		plot.title = element_text(size=10),\n";
	rSlopeBoxPlt << "		axis.text = element_text(size = 8, colour = \"black\"),\n";
	rSlopeBoxPlt << "		axis.text.x=element_blank(),\n";
	rSlopeBoxPlt << "		legend.key = element_rect(fill = \"white\"),\n";
	rSlopeBoxPlt << "		legend.background = element_rect(fill = \"white\"),\n";
	rSlopeBoxPlt << "		legend.position = (\"bottom\"),\n";
	rSlopeBoxPlt << "		legend.text = element_text(size = 7),\n";
	rSlopeBoxPlt << "		legend.title=element_blank(),\n";
	rSlopeBoxPlt << "		panel.grid.major = element_line(colour = \"white\"),\n";
	rSlopeBoxPlt << "		panel.grid.minor = element_blank(),\n";
	rSlopeBoxPlt << "		panel.background = element_rect(fill = \"white\"),\n";
	rSlopeBoxPlt << "		axis.line.x = element_line(colour = \"black\"),\n";
	rSlopeBoxPlt << "		axis.line.y = element_line(colour = \"black\")\n";
	rSlopeBoxPlt << "	)\n";


	rSlopeBoxPlt.close();
	csvSlope.close();
}




