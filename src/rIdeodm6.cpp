#include <iostream>
#include <fstream>
#include <string>
//#include <regex> 
#include <vector>
#include <algorithm>
using namespace std;

std::string exec(const char* cmd);

// creates ideogram with trunover heatmap
void rIdeodm6(std::string s_rPltsName, int i_bamFiles, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, std::string s_direction)  
{
	// write r script
	ofstream rIdeo;
	rIdeo.open (s_rPltsName, std::ios::app);

	// write csv files // XXX ADDED XXX
	ofstream csvX1;
	ofstream csvX2;
	ofstream csvY1;
	ofstream csvY2;
	ofstream csvInf;
	std::string s_csvX1Name;
	std::string s_csvX2Name;
	std::string s_csvY1Name;
	std::string s_csvY2Name;
	std::string s_csvInfName;
	if (s_direction == "forward") { s_csvX1Name = "./" + s_rPltsName + "_data/fIdeoX1.csv"; }
	if (s_direction == "reverse") { s_csvX1Name = "./" + s_rPltsName + "_data/rIdeoX1.csv"; }
	if (s_direction == "forward") { s_csvX2Name = "./" + s_rPltsName + "_data/fIdeoX2.csv"; }
	if (s_direction == "reverse") { s_csvX2Name = "./" + s_rPltsName + "_data/rIdeoX2.csv"; }
	if (s_direction == "forward") { s_csvY1Name = "./" + s_rPltsName + "_data/fIdeoY1.csv"; }
	if (s_direction == "reverse") { s_csvY1Name = "./" + s_rPltsName + "_data/rIdeoY1.csv"; }
	if (s_direction == "forward") { s_csvY2Name = "./" + s_rPltsName + "_data/fIdeoY2.csv"; }
	if (s_direction == "reverse") { s_csvY2Name = "./" + s_rPltsName + "_data/rIdeoY2.csv"; }
	if (s_direction == "forward") { s_csvInfName = "./" + s_rPltsName + "_data/fIdeoInf.csv"; }
	if (s_direction == "reverse") { s_csvInfName = "./" + s_rPltsName + "_data/rIdeoInf.csv"; }
	csvX1.open(s_csvX1Name);
	csvX2.open(s_csvX2Name);
	csvY1.open(s_csvY1Name);
	csvY2.open(s_csvY2Name);
	csvInf.open(s_csvInfName);

	int i_goodpeak = 0;

	// count non NAN containing turnover cells and set count to i_goodpeak
	for (int i = 1; i < i_peakNumber + 1; i++) { 
		if (s_direction == "forward") {					
			std::string s = dataArray[i][2*i_bamFiles + 4];		
			if ( (s != "NaN") && (s != "-") ) { i_goodpeak++; }	
		}
		if (s_direction == "reverse") {					
			std::string s = dataArray[i][2*i_bamFiles + 8];		
			if ( (s != "NaN") && (s != "-") ) { i_goodpeak++; }	
		}				
	}


	// create array of non NAN containing turnover chromosomes and their start and end sites
	std::string chrArr[i_goodpeak];
	int chrStart[i_goodpeak];
	int chrEnd[i_goodpeak];

	int i_iterator = 0;
	for (int i = 1; i < i_peakNumber + 1; i++) {
		if (s_direction == "forward") {							
			std::string s_turnover = dataArray[i][2*i_bamFiles + 4];		
			if ( (s_turnover != "NaN") && (s_turnover != "-") ) {			
				std::string s_chrCell = dataArray[i][0];		
				s_chrCell.erase(0, 3);
				chrArr[i_iterator] = s_chrCell;
				chrStart[i_iterator] = stoi(dataArray[i][1].c_str());
				chrEnd[i_iterator] = stoi(dataArray[i][2].c_str());
				i_iterator++;
			}
		}
		if (s_direction == "reverse") {							
			std::string s_turnover = dataArray[i][2*i_bamFiles + 8];		
			if ( (s_turnover != "NaN") && (s_turnover != "-") ) {			
				std::string s_chrCell = dataArray[i][0];		
				s_chrCell.erase(0, 3);
				chrArr[i_iterator] = s_chrCell;
				chrStart[i_iterator] = stoi(dataArray[i][1].c_str());
				chrEnd[i_iterator] = stoi(dataArray[i][2].c_str());
				i_iterator++;
			}
		}
	}



	rIdeo << "#dm6\n";
	if (s_direction == "forward") {							
		rIdeo << "fIdeoX1 <- suppressWarnings(read.csv(file=\"./" + s_rPltsName + "_data/fIdeoX1.csv\", header=FALSE))\n";
		rIdeo << "fIdeoX2 <- suppressWarnings(read.csv(file=\"./" + s_rPltsName + "_data/fIdeoX2.csv\", header=FALSE))\n";
		rIdeo << "fIdeoY1 <- suppressWarnings(read.csv(file=\"./" + s_rPltsName + "_data/fIdeoY1.csv\", header=FALSE))\n";
		rIdeo << "fIdeoY2 <- suppressWarnings(read.csv(file=\"./" + s_rPltsName + "_data/fIdeoY2.csv\", header=FALSE))\n";
		rIdeo << "fIdeoInf <- suppressWarnings(read.csv(file=\"./" + s_rPltsName + "_data/fIdeoInf.csv\", header=FALSE))\n";
		rIdeo << "fIdeo <- data.frame(fIdeoX1, fIdeoX2, fIdeoY1, fIdeoY2, fIdeoInf)\n";
	}
	if (s_direction == "reverse") {	
		rIdeo << "rIdeoX1 <- suppressWarnings(read.csv(file=\"./" + s_rPltsName + "_data/rIdeoX1.csv\", header=FALSE))\n";
		rIdeo << "rIdeoX2 <- suppressWarnings(read.csv(file=\"./" + s_rPltsName + "_data/rIdeoX2.csv\", header=FALSE))\n";
		rIdeo << "rIdeoY1 <- suppressWarnings(read.csv(file=\"./" + s_rPltsName + "_data/rIdeoY1.csv\", header=FALSE))\n";
		rIdeo << "rIdeoY2 <- suppressWarnings(read.csv(file=\"./" + s_rPltsName + "_data/rIdeoY2.csv\", header=FALSE))\n";
		rIdeo << "rIdeoInf <- suppressWarnings(read.csv(file=\"./" + s_rPltsName + "_data/rIdeoInf.csv\", header=FALSE))\n";
		rIdeo << "rIdeo <- data.frame(rIdeoX1, rIdeoX2, rIdeoY1, rIdeoY2, rIdeoInf)\n";
	}

	// REQUIRED: chromosome, start, end, turnover, # of valid turnover peaks.
	// VARIABLE PORTION XXX

	// ideox1 position of band
	for (int i = 0; i < i_goodpeak-1; i++) {
		if (chrArr[i] == "3R") {
			csvX1 << "0\n";
			csvX2 << "1\n";
		} else if (chrArr[i] == "3L") {
			csvX1 << "1.5\n";
			csvX2 << "2.5\n";
		} else if (chrArr[i] == "2R") {
			csvX1 << "3\n";
			csvX2 << "4\n";
		} else if (chrArr[i] == "X") {
			csvX1 << "4.5\n";
			csvX2 << "5.5\n";
		} else if (chrArr[i] == "2L") {
			csvX1 << "6\n";
			csvX2 << "7\n";
		} else if (chrArr[i] == "Y") {
			csvX1 << "7.5\n";
			csvX2 << "8.5\n";
		} else if (chrArr[i] == "4") {
			csvX1 << "9\n";
			csvX2 << "10\n";
		} else { // different chromosomes written off grid boundaries
			csvX1 << "-5\n";
			csvX2 << "-5\n";
		} 
	}

	// ideox1 end position of band
	if (chrArr[i_goodpeak-1] == "3R") {
		csvX1 << "0";
		csvX2 << "1";
	} else if (chrArr[i_goodpeak - 1] == "3L") {
		csvX1 << "1.5";
		csvX2 << "2.5";
	} else if (chrArr[i_goodpeak - 1] == "2R") {
		csvX1 << "3";
		csvX2 << "4";
	} else if (chrArr[i_goodpeak - 1] == "X") {
		csvX1 << "4.5";
		csvX2 << "5.5";
	} else if (chrArr[i_goodpeak - 1] == "2L") {
		csvX1 << "6";
		csvX2 << "7";
	} else if (chrArr[i_goodpeak - 1] == "Y") {
		csvX1 << "7.5";
		csvX2 << "8.5";
	} else if (chrArr[i_goodpeak - 1] == "4") {
		csvX1 << "9";
		csvX2 << "10";
	} else { // different chromosomes written off grid boundaries
		csvX1 << "-5";
		csvX2 << "-5";
	}


	// band normalized from center to 10000 bp XXX

	// ideoy1 position of band
	int i_mid = 0;
	for (int i = 0; i < i_goodpeak-1; i++) {
		i_mid = (chrEnd[i] - chrStart[i])/2;
		csvY1 << chrStart[i] + i_mid - 10000 << "\n";
	}
	i_mid = (chrEnd[i_goodpeak] - chrStart[i_goodpeak])/2;
	csvY1 << chrStart[i_goodpeak] + i_mid - 10000;

	// ideoy2 position of band
	i_mid = 0;
	for (int i = 0; i < i_goodpeak-1; i++) {
		i_mid = (chrEnd[i] - chrStart[i])/2;
		csvY2 << chrStart[i] + i_mid + 10000 << "\n";
	}
	i_mid = (chrEnd[i_goodpeak] - chrStart[i_goodpeak])/2;
	csvY2 << chrStart[i_goodpeak] + i_mid + 10000;

	i_iterator = 0;
	// t values of bands
	for (int i = 1; i < i_peakNumber + 1; i++) { 							
		if (s_direction == "forward") {								
			std::string s = dataArray[i][2*i_bamFiles + 4];					
			if ( (s != "NaN") && (s != "-") && (i_iterator == i_goodpeak-1) ) {	
				csvInf << s;								
			}												
			if ( (s != "NaN") && (s != "-") && (i_iterator != i_goodpeak-1) ) {	
				csvInf << s << "\n";								
				i_iterator++;									
			}												
		}													
		if (s_direction == "reverse") {								
			std::string s = dataArray[i][2*i_bamFiles + 8];					
			if ( (s != "NaN") && (s != "-") && (i_iterator == i_goodpeak-1) ) {	
				csvInf << s;								
			}	
			if ( (s != "NaN") && (s != "-") && (i_iterator != i_goodpeak-1) ) {	
				csvInf << s << "\n";								
				i_iterator++;									
			}													
		}
	}

	// Draw all chromosomes

	rIdeo << "ch=data.frame(chx=c(0.5,2,3.5,5,6.5,8,9.5), t=c(\"3R\",\"3L\",\"2R\",\"X\",\"2L\",\"Y\",\"4\"))\n";
	rIdeo << "d=data.frame(x=c(0,1,1,0, 1.5,2.5,2.5,1.5, 3,4,4,3, 4.5,5.5,5.5,4.5, 6,7,7,6, 7.5,8.5,8.5,7.5, 9,10,10,9), y=c(0,0,32079331,32079331, 0,0,28110227,28110227, 0,0,25286936,25286936, 0,0,23542271,23542271, 0,0,23513712,23513712, 0,0,3667352,3667352, 0,0,1348131,1348131), t=c('1', '1', '1', '1', '2', '2', '2', '2', '3', '3', '3', '3', '4', '4', '4', '4', '5', '5', '5', '5', '6', '6', '6', '6', '7', '7', '7', '7'))\n";

	// this is same for all except the scale_x_continuous(limits = c(0, 32)) - 32 is variable

	rIdeo << "#y=c(cenS,cenE,chrE,chrE,cenE,cenS,0,0,\n";
	if (s_direction == "forward") { rIdeo << "ideForward <- ggplot() +\n"; }			
	if (s_direction == "reverse") { rIdeo << "ideReverse <- ggplot() +\n"; }			
	rIdeo << "  scale_x_continuous(limits = c(0, 10.5)) +\n";
	rIdeo << "  scale_y_continuous(limits = c(-5000000, 33000000)) +\n";
	rIdeo << "  geom_polygon(data=d, mapping=aes(x=x, y=y, group=t), colour=\"black\", fill=\"white\") + \n";
	if (s_direction == "forward") {
  		rIdeo << "geom_rect(data=fIdeo, mapping=aes(xmin=fIdeo[,1], xmax=fIdeo[,2], ymin=fIdeo[,3], ymax=fIdeo[,4], fill=fIdeo[,5])) +\n";
	} else {
  		rIdeo << "geom_rect(data=rIdeo, mapping=aes(xmin=rIdeo[,1], xmax=rIdeo[,2], ymin=rIdeo[,3], ymax=rIdeo[,4], fill=rIdeo[,5])) +\n";
	}
	rIdeo << "  geom_polygon(data=d, mapping=aes(x=x, y=y, group=t), colour=\"black\", alpha=0) +\n";
	rIdeo << "  scale_fill_gradientn(trans = \"log10\",colours=rainbow(14)) +\n";
	rIdeo << "  geom_text(data=ch, aes(x=chx, y=-5000000, label=t), size=4) +\n";
	if (s_direction == "forward") { 
		rIdeo << "  ggtitle(\"Ideogram Heat Map of TTI at Signal Increase Loci\") +\n"; }	
	if (s_direction == "reverse") { 
		rIdeo << "  ggtitle(\"Ideogram Heat Map of TTI at Signal Decrease Loci\") +\n"; }	
	rIdeo << "  labs(x = \"Chromosome\",size = 14) +\n";
	rIdeo << "  labs(y = \"Size (bp)\",size = 14) +\n";
	rIdeo << "  labs(fill = \"\") +\n";
	rIdeo << "  theme(\n";
	rIdeo << "    axis.text = element_text(size = 8, colour = \"black\"),\n";
	rIdeo << "    legend.key = element_rect(fill = \"white\"),\n";
	rIdeo << "    legend.background = element_rect(fill = \"white\"),\n";
	rIdeo << "    panel.grid.major = element_line(colour = \"white\"),\n";
	rIdeo << "    panel.grid.minor = element_blank(),\n";
	rIdeo << "    panel.background = element_rect(fill = \"grey95\"),\n";
	rIdeo << "    legend.key.height=unit(2,\"cm\"),\n";
	rIdeo << "    axis.text.x = element_blank(),\n";
	rIdeo << "    axis.ticks.x = element_blank(),\n";
	rIdeo << "    #axis.line.x = element_line(colour = \"black\"),\n";
	rIdeo << "    axis.line.y = element_line(colour = \"black\")\n";
	rIdeo << "  )\n";

	rIdeo.close();
	csvX1.close();
	csvX2.close();
	csvY1.close();
	csvY2.close();
	csvInf.close();
}

