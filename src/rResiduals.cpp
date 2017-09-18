#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;


void rResiduals(std::string s_rPltsName, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, int i_bamFiles) 
{

	// write r script
	ofstream rResBoxPlt;
	rResBoxPlt.open(s_rPltsName, std::ios::app);

	// write csv file
	ofstream csvRes;
	ofstream csvID;
	std::string s_ResCsvName;
	std::string s_idCsvName;
	
	s_ResCsvName = "./" + s_rPltsName + "_data/ResBox.csv";
	s_idCsvName = "./" + s_rPltsName + "_data/idRes.csv";

	csvRes.open(s_ResCsvName);
	csvID.open(s_idCsvName);

	rResBoxPlt << "ResBox <- read.csv(file=\"./" << s_rPltsName << "_data/ResBox.csv\", header=FALSE)\n";
	rResBoxPlt << "idRes <- read.csv(file=\"./" << s_rPltsName << "_data/idRes.csv\", header=FALSE, stringsAsFactors=FALSE)\n";


	for (int j = 1; j < i_peakNumber; j++) { 
		std::string s_type = dataArray[j][2*i_bamFiles + 4]; // locus type
		std::string s_riseRes =  dataArray[j][2*i_bamFiles + 11] + '\n';
		std::string s_fallRes =  dataArray[j][2*i_bamFiles + 18] + '\n';
		if (s_type == "rise") {  csvID << s_type << '\n'; csvRes << s_riseRes;}
		if (s_type == "fall") {  csvID << s_type << '\n'; csvRes << s_fallRes;}
		if (s_type == "undefined") {  
			if (s_riseRes == "-") { // undefined fall
				csvID << "undefined fall\n"; csvRes << s_fallRes;
			} else { // undefined rise
				csvID << "undefined rise\n"; csvRes << s_riseRes;
			}
		}
		if (s_type == "hill") { csvID << "hill incline\nhill decline\n"; csvRes << s_riseRes << s_fallRes; }
		if (s_type == "valley") { csvID << "valley decline\nvalley incline\n"; csvRes << s_fallRes << s_riseRes; }
	}

	std::string s_type = dataArray[i_peakNumber][2*i_bamFiles + 4]; // locus type
	if (s_type == "rise") {  csvID << s_type; csvRes << dataArray[i_peakNumber][2*i_bamFiles + 11];}
	if (s_type == "fall") {  csvID << s_type; csvRes << dataArray[i_peakNumber][2*i_bamFiles + 18];}
	if (s_type == "undefined") {  
		if (dataArray[i_peakNumber][2*i_bamFiles + 11] == "-") { // undefined fall
			csvID << "undefined fall"; csvRes << dataArray[i_peakNumber][2*i_bamFiles + 18];
		} else { // undefined rise
			csvID << "undefined rise"; csvRes << dataArray[i_peakNumber][2*i_bamFiles + 11];
		}
	}
	if (s_type == "hill") { csvID << "hill incline\nhill decline"; csvRes << dataArray[i_peakNumber][2*i_bamFiles + 11] << "\n" << dataArray[i_peakNumber][2*i_bamFiles + 18]; }
	if (s_type == "valley") { csvID << "valley decline\nvalley incline"; csvRes << dataArray[i_peakNumber][2*i_bamFiles + 18] << "\n" << dataArray[i_peakNumber][2*i_bamFiles + 11]; }

	
	rResBoxPlt << "dfResBox<- data.frame(idRes,ResBox)\n"; 
	rResBoxPlt << "sigRes <- ggplot(dfResBox, aes(dfResBox[,1],dfResBox[,2], fill=dfResBox[,1])) +\n"; 
	rResBoxPlt << "	geom_boxplot(outlier.shape = NA) +\n";
	rResBoxPlt << "	scale_fill_brewer(palette=\"Dark2\") +\n";
	rResBoxPlt << "	ggtitle(\"Residuals from Sigmoidal Fits by Locus Type\") +\n";


	rResBoxPlt << "	labs(y = \"Sigmoidal Residuals\",size = 14) +\n";
	rResBoxPlt << "	labs(x = \"Locus Type\",size = 14) +\n";
	rResBoxPlt << "	theme(\n";
	rResBoxPlt << "		axis.text = element_text(size = 8, colour = \"black\"),\n";
	rResBoxPlt << "		plot.title = element_text(size=10),\n";
	rResBoxPlt << "		axis.text.x=element_blank(),\n";
	rResBoxPlt << "		legend.key = element_rect(fill = \"white\"),\n";
	rResBoxPlt << "		legend.background = element_rect(fill = \"white\"),\n";
	rResBoxPlt << "		legend.position = (\"bottom\"),\n";
	rResBoxPlt << "		legend.text = element_text(size = 7),\n";
	rResBoxPlt << "		legend.title=element_blank(),\n";
	rResBoxPlt << "		panel.grid.major = element_line(colour = \"white\"),\n";
	rResBoxPlt << "		panel.grid.minor = element_blank(),\n";
	rResBoxPlt << "		panel.background = element_rect(fill = \"white\"),\n";
	rResBoxPlt << "		axis.line.x = element_line(colour = \"black\"),\n";
	rResBoxPlt << "		axis.line.y = element_line(colour = \"black\")\n";
	rResBoxPlt << "	)\n";


	rResBoxPlt.close();
	csvRes.close();
	csvID.close();
}




