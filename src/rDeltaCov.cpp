#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;


void rDeltaCov(std::string s_rPltsName, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, int i_bamFiles) 
{

	// write r script
	ofstream rDCovBoxPlt;
	rDCovBoxPlt.open(s_rPltsName, std::ios::app);

	// write csv file
	ofstream csvCov;
	ofstream csvID;
	std::string s_dCovCsvName;
	std::string s_idCsvName;
	
	s_dCovCsvName = "./" + s_rPltsName + "_data/dCovBox.csv";
	s_idCsvName = "./" + s_rPltsName + "_data/idDCov.csv";

	csvCov.open(s_dCovCsvName);
	csvID.open(s_idCsvName);

	rDCovBoxPlt << "dCovBox <- read.csv(file=\"./" << s_rPltsName << "_data/dCovBox.csv\", header=FALSE)\n";
	rDCovBoxPlt << "idDCov <- read.csv(file=\"./" << s_rPltsName << "_data/idDCov.csv\", header=FALSE, stringsAsFactors=FALSE)\n";

	for (int j = 1; j < i_peakNumber; j++) { 
		std::string s =  dataArray[j][2*i_bamFiles + 3] + '\n';
		csvCov << s;  // delta coverage
		s = dataArray[j][2*i_bamFiles + 4] + '\n'; // locus type
		csvID << s;
	}
	csvCov << dataArray[i_peakNumber][2*i_bamFiles + 3];
	csvID << dataArray[i_peakNumber][2*i_bamFiles + 4];


	
	rDCovBoxPlt << "dfCovBox<- data.frame(idDCov,dCovBox)\n"; 
	rDCovBoxPlt << "dCov <- ggplot(dfCovBox, aes(dfCovBox[,1],dfCovBox[,2], fill=dfCovBox[,1])) +\n"; 
	rDCovBoxPlt << "	geom_boxplot(outlier.shape = NA) +\n";
	rDCovBoxPlt << "	scale_fill_brewer(palette=\"Dark2\") +\n";
	rDCovBoxPlt << "	ggtitle(\"Delta Coverage by Locus Type\") +\n";


	rDCovBoxPlt << "	labs(y = \"Delta Coverage\",size = 14) +\n";
	rDCovBoxPlt << "	labs(x = \"Locus Type\",size = 14) +\n";
	rDCovBoxPlt << "	theme(\n";
	rDCovBoxPlt << "		plot.title = element_text(size=10),\n";
	rDCovBoxPlt << "		axis.text = element_text(size = 8, colour = \"black\"),\n";
	rDCovBoxPlt << "		axis.text.x=element_blank(),\n";
	rDCovBoxPlt << "		legend.key = element_rect(fill = \"white\"),\n";
	rDCovBoxPlt << "		legend.background = element_rect(fill = \"white\"),\n";
	rDCovBoxPlt << "		legend.position = (\"bottom\"),\n";
	rDCovBoxPlt << "		legend.text = element_text(size = 7),\n";
	rDCovBoxPlt << "		legend.title=element_blank(),\n";
	rDCovBoxPlt << "		panel.grid.major = element_line(colour = \"white\"),\n";
	rDCovBoxPlt << "		panel.grid.minor = element_blank(),\n";
	rDCovBoxPlt << "		panel.background = element_rect(fill = \"white\"),\n";
	rDCovBoxPlt << "		axis.line.x = element_line(colour = \"black\"),\n";
	rDCovBoxPlt << "		axis.line.y = element_line(colour = \"black\")\n";
	rDCovBoxPlt << "	)\n";


	rDCovBoxPlt.close();
	csvCov.close();
	csvID.close();
}




