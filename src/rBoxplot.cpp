#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;


void rBoxplot(int turnoverTimes[], std::string s_rPltsName, int i_peakNumber, int i_genomeFileCount, int i_bamFiles, std::vector<std::vector<std::string> > &dataArray, std::string s_direction) 
{

	// write r script
	ofstream rBoxplotPlt;
	rBoxplotPlt.open(s_rPltsName, std::ios::app);

	// write csv file
	ofstream csvInf;
	ofstream csvID;
	std::string s_infCsvName;
	std::string s_idCsvName;
	if (s_direction == "forward") { s_infCsvName = "./" + s_rPltsName + "_data/fInfBox.csv"; }
	if (s_direction == "reverse") { s_infCsvName = "./" + s_rPltsName + "_data/rInfBox.csv"; }
	if (s_direction == "forward") { s_idCsvName = "./" + s_rPltsName + "_data/idBoxForward.csv"; }
	if (s_direction == "reverse") { s_idCsvName = "./" + s_rPltsName + "_data/idBoxReverse.csv"; }
	csvInf.open(s_infCsvName);
	csvID.open(s_idCsvName);


	if (s_direction == "forward") { 
		rBoxplotPlt << "fInfBox <- read.csv(file=\"./" << s_rPltsName << "_data/fInfBox.csv\", header=FALSE)\n";
	} else { // reverse
		rBoxplotPlt << "rInfBox <- read.csv(file=\"./" << s_rPltsName << "_data/rInfBox.csv\", header=FALSE)\n";
	}

	int i_arr[i_genomeFileCount]; // valid overlaps per gene feature
	int i_colstart = 2*i_bamFiles + 19;
	int i_colend = 2*i_bamFiles + 19 + i_genomeFileCount;
	int i_loop = 0;

	for (int i = i_colstart; i < i_colend; i++) { // for each gene feature
		int i_iterator = 0;
		for (int j = 1; j < i_peakNumber + 1; j++) { // number of valid turnover values for each gene feature
			if (s_direction == "forward") {					
				std::string s = dataArray[j][2*i_bamFiles + 5];	// rise TTI
				if ( (atoi(dataArray[j][i].c_str()) > 0) && (s != "NaN") && (s != "-") ) {
					i_iterator++;
				}
			}
			if (s_direction == "reverse") {					
				std::string s = dataArray[j][2*i_bamFiles + 12];	// fall TTI
				if ( (atoi(dataArray[j][i].c_str()) > 0) && (s != "NaN") && (s != "-") ) {
					i_iterator++;
				}
			}
		}
		i_arr[i_loop] = i_iterator;
		i_loop++;
	}

	// inflection points
	for (int i = i_colstart; i < i_colend-1; i++) { // for each gene feature except last
		for (int j = 1; j < i_peakNumber + 1; j++) {
			std::string s;
			if (s_direction == "forward") { s = dataArray[j][2*i_bamFiles + 5]; } // rise TTI				
			if (s_direction == "reverse") { s = dataArray[j][2*i_bamFiles + 12]; } // fall TTI				
			if ( (atoi(dataArray[j][i].c_str()) > 0) && (s != "NaN") && (s != "-") ) { csvInf << s << "\n"; }
		}
	}

	int i_iterator = 1;
	for (int j = 1; j < i_peakNumber + 1; j++) { // for turnover values of last gene feature	
		std::string s;
		if (s_direction == "forward") { s = dataArray[j][2*i_bamFiles + 5]; } // rise TTI
		if (s_direction == "reverse") { s = dataArray[j][2*i_bamFiles + 12]; } // fall TTI
		// put NA value if last genome feature has zero overlaps
		if (i_arr[i_genomeFileCount-1] == 0) {
			csvInf << "NA";
			break;
		}
		if ( (atoi(dataArray[j][i_colend-1].c_str()) > 0) && (s != "NaN") && (s != "-") && (i_arr[i_genomeFileCount-1]) > (i_iterator) ) {
			csvInf << s << "\n";
			i_iterator++;
		} else if ( (atoi(dataArray[j][i_colend-1].c_str()) > 0) && (s != "NaN") && (s != "-") && (i_arr[i_genomeFileCount-1] == (i_iterator)) ) {
			csvInf << s;
		} else {}
	}


	// ID
	if (s_direction == "forward") { 
		rBoxplotPlt << "fIDBox <- read.csv(file=\"./" << s_rPltsName << "_data/idBoxForward.csv\", header=FALSE, stringsAsFactors=FALSE)\n";
	} else { // reverse
		rBoxplotPlt << "rIDBox <- read.csv(file=\"./" << s_rPltsName << "_data/idBoxReverse.csv\", header=FALSE, stringsAsFactors=FALSE)\n";
	}

	for (int i = 0; i < i_genomeFileCount-1; i++) {
		for (int j = 0; j < i_arr[i]; j++) {
			std::string s_feat = dataArray[0][i_colstart+i];
			s_feat = s_feat.substr(0, s_feat.size()-4);
			csvID << s_feat << "\n";
		}
	}

	for (int i = 0; i < i_arr[i_genomeFileCount-1]-1; i++) {
		std::string s_feat = dataArray[0][i_colend-1];
		s_feat = s_feat.substr(0, s_feat.size()-4);
		csvID << s_feat << "\n";
	}
	std::string s_feat = dataArray[0][i_colend-1];
	s_feat = s_feat.substr(0, s_feat.size()-4);
	csvID << s_feat;

	
	if (s_direction == "forward") { rBoxplotPlt << "fBox<- data.frame(fIDBox,fInfBox)\n"; }
	if (s_direction == "reverse") { rBoxplotPlt << "rBox<- data.frame(rIDBox,rInfBox)\n"; }
	if (s_direction == "forward") { rBoxplotPlt << "boxForward <- ggplot(fBox, aes(fBox[,1],fBox[,2], fill=fBox[,1])) +\n"; }
	if (s_direction == "reverse") { rBoxplotPlt << "boxReverse <- ggplot(rBox, aes(rBox[,1],rBox[,2], fill=rBox[,1])) +\n"; }
	rBoxplotPlt << "	geom_boxplot(outlier.shape = NA) +\n";
	rBoxplotPlt << "	scale_fill_brewer(palette=\"Dark2\") +\n";
	if (s_direction == "forward") {
	rBoxplotPlt << "	ggtitle(\"TTI of Signal Increase Loci at Gene Features\") +\n"; }
	if (s_direction == "reverse") {
	rBoxplotPlt << "	ggtitle(\"TTI of Signal Decrease Loci at Gene Features\") +\n"; }
	rBoxplotPlt << "	labs(y = \"TTI (time)\",size = 14) +\n";
	rBoxplotPlt << "	labs(x = \"Feature\",size = 14) +\n";

	rBoxplotPlt << "	scale_y_continuous(limits = c(0, " << turnoverTimes[i_bamFiles-1]*2 << ")) +\n";
	rBoxplotPlt << "	theme(\n";
	rBoxplotPlt << "		axis.text = element_text(size = 8, colour = \"black\"),\n";
	rBoxplotPlt << "		axis.text.x=element_blank(),\n";
	rBoxplotPlt << "		legend.key = element_rect(fill = \"white\"),\n";
	rBoxplotPlt << "		legend.background = element_rect(fill = \"white\"),\n";
	rBoxplotPlt << "		legend.position = (\"bottom\"),\n";
	rBoxplotPlt << "		legend.text = element_text(size = 7),\n";
	rBoxplotPlt << "		legend.title=element_blank(),\n";
	rBoxplotPlt << "		panel.grid.major = element_line(colour = \"white\"),\n";
	rBoxplotPlt << "		panel.grid.minor = element_blank(),\n";
	rBoxplotPlt << "		panel.background = element_rect(fill = \"white\"),\n";
	rBoxplotPlt << "		axis.line.x = element_line(colour = \"black\"),\n";
	rBoxplotPlt << "		axis.line.y = element_line(colour = \"black\")\n";
	rBoxplotPlt << "	)\n";


	rBoxplotPlt.close();
	csvInf.close();
	csvID.close();
}




