#include <iostream>
#include <fstream>
#include <string>
//#include <regex> 
#include <vector>
#include <algorithm>
using namespace std;

std::string exec(const char* cmd);

// creates ideogram with trunover heatmap
void rIdeomm10(std::string s_rPltsName, int i_bamFiles, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, std::string s_direction)  
{
	// write r script
	ofstream rIdeo;
	rIdeo.open (s_rPltsName, std::ios::app);

	// write csv files
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



	rIdeo << "#mm10\n";
	if (s_direction == "forward") {							
		rIdeo << "fIdeoX1 <- read.csv(file=\"./" + s_rPltsName + "_data/fIdeoX1.csv\", header=FALSE)\n";
		rIdeo << "fIdeoX2 <- read.csv(file=\"./" + s_rPltsName + "_data/fIdeoX2.csv\", header=FALSE)\n";
		rIdeo << "fIdeoY1 <- read.csv(file=\"./" + s_rPltsName + "_data/fIdeoY1.csv\", header=FALSE)\n";
		rIdeo << "fIdeoY2 <- read.csv(file=\"./" + s_rPltsName + "_data/fIdeoY2.csv\", header=FALSE)\n";
		rIdeo << "fIdeoInf <- read.csv(file=\"./" + s_rPltsName + "_data/fIdeoInf.csv\", header=FALSE)\n";
		rIdeo << "fIdeo <- data.frame(fIdeoX1, fIdeoX2, fIdeoY1, fIdeoY2, fIdeoInf)\n";
	}
	if (s_direction == "reverse") {	
		rIdeo << "rIdeoX1 <- read.csv(file=\"./" + s_rPltsName + "_data/rIdeoX1.csv\", header=FALSE)\n";
		rIdeo << "rIdeoX2 <- read.csv(file=\"./" + s_rPltsName + "_data/rIdeoX2.csv\", header=FALSE)\n";
		rIdeo << "rIdeoY1 <- read.csv(file=\"./" + s_rPltsName + "_data/rIdeoY1.csv\", header=FALSE)\n";
		rIdeo << "rIdeoY2 <- read.csv(file=\"./" + s_rPltsName + "_data/rIdeoY2.csv\", header=FALSE)\n";
		rIdeo << "rIdeoInf <- read.csv(file=\"./" + s_rPltsName + "_data/rIdeoInf.csv\", header=FALSE)\n";
		rIdeo << "rIdeo <- data.frame(rIdeoX1, rIdeoX2, rIdeoY1, rIdeoY2, rIdeoInf)\n";
	}

	// REQUIRED: chromosome, start, end, turnover, # of valid turnover peaks.
	// VARIABLE PORTION XXX
	// x position of band
	for (int i = 0; i < i_goodpeak-1; i++) {
		if (chrArr[i] == "1") {
			csvX1 << "0\n";
			csvX2 << "1\n";
		} else if (chrArr[i] == "2") {
			csvX1 << "1.5\n";
			csvX2 << "2.5\n";
		} else if (chrArr[i] == "3") {
			csvX1 << "3\n";
			csvX2 << "4\n";
		} else if (chrArr[i] == "4") {
			csvX1 << "4.5\n";
			csvX2 << "5.5\n";
		} else if (chrArr[i] == "5") {
			csvX1 << "6\n";
			csvX2 << "7\n";
		} else if (chrArr[i] == "6") {
			csvX1 << "7.5\n";
			csvX2 << "8.5\n";
		} else if (chrArr[i] == "7") {
			csvX1 << "9\n";
			csvX2 << "10\n";
		} else if (chrArr[i] == "8") {
			csvX1 << "10.5\n";
			csvX2 << "11.5\n";
		} else if (chrArr[i] == "9") {
			csvX1 << "12\n";
			csvX2 << "13\n";
		} else if (chrArr[i] == "10") {
			csvX1 << "13.5\n";
			csvX2 << "14.5\n";
		} else if (chrArr[i] == "11") {
			csvX1 << "15\n";
			csvX2 << "16\n";
		} else if (chrArr[i] == "12") {
			csvX1 << "16.5\n";
			csvX2 << "17.5\n";
		} else if (chrArr[i] == "13") {
			csvX1 << "18\n";
			csvX2 << "19\n";
		} else if (chrArr[i] == "14") {
			csvX1 << "19.5\n";
			csvX2 << "20.5\n";
		} else if (chrArr[i] == "15") {
			csvX1 << "21\n";
			csvX2 << "22\n";
		} else if (chrArr[i] == "16") {
			csvX1 << "22.5\n";
			csvX2 << "23.5\n";
		} else if (chrArr[i] == "17") {
			csvX1 << "24\n";
			csvX2 << "25\n";
		} else if (chrArr[i] == "18") {
			csvX1 << "25.5\n";
			csvX2 << "26.5\n";
		} else if (chrArr[i] == "19") {
			csvX1 << "27\n";
			csvX2 << "28\n";
		} else if (chrArr[i] == "X") {
			csvX1 << "28.5\n";
			csvX2 << "29.5\n";
		} else if (chrArr[i] == "Y") {
			csvX1 << "30\n";
			csvX2 << "31\n";
		} else { // different chromosomes written off grid boundaries
			csvX1 << "-5\n";
			csvX2 << "-5\n";
		} 
	}

	// last x position of band
	if (chrArr[i_goodpeak-1] == "1") {
		csvX1 << "0";
		csvX2 << "1";
	} else if (chrArr[i_goodpeak - 1] == "2") {
		csvX1 << "1.5";
		csvX2 << "2.5";
	} else if (chrArr[i_goodpeak - 1] == "3") {
		csvX1 << "3";
		csvX2 << "4";
	} else if (chrArr[i_goodpeak - 1] == "4") {
		csvX1 << "4.5";
		csvX2 << "5.5";
	} else if (chrArr[i_goodpeak - 1] == "5") {
		csvX1 << "6";
		csvX2 << "7";
	} else if (chrArr[i_goodpeak - 1] == "6") {
		csvX1 << "7.5";
		csvX2 << "8.5";
	} else if (chrArr[i_goodpeak - 1] == "7") {
		csvX1 << "9";
		csvX2 << "10";
	} else if (chrArr[i_goodpeak - 1] == "8") {
		csvX1 << "10.5";
		csvX2 << "11.5";
	} else if (chrArr[i_goodpeak - 1] == "9") {
		csvX1 << "12";
		csvX2 << "13";
	} else if (chrArr[i_goodpeak - 1] == "10") {
		csvX1 << "13.5";
		csvX2 << "14.5";
	} else if (chrArr[i_goodpeak - 1] == "11") {
		csvX1 << "15";
		csvX2 << "16";
	} else if (chrArr[i_goodpeak - 1] == "12") {
		csvX1 << "16.5";
		csvX2 << "17.5";
	} else if (chrArr[i_goodpeak - 1] == "13") {
		csvX1 << "18";
		csvX2 << "19";
	} else if (chrArr[i_goodpeak - 1] == "14") {
		csvX1 << "19.5";
		csvX2 << "20.5";
	} else if (chrArr[i_goodpeak - 1] == "15") {
		csvX1 << "21";
		csvX2 << "22";
	} else if (chrArr[i_goodpeak - 1] == "16") {
		csvX1 << "22.5";
		csvX2 << "23.5";
	} else if (chrArr[i_goodpeak - 1] == "17") {
		csvX1 << "24";
		csvX2 << "25";
	} else if (chrArr[i_goodpeak - 1] == "18") {
		csvX1 << "25.5";
		csvX2 << "26.5";
	} else if (chrArr[i_goodpeak - 1] == "19") {
		csvX1 << "27";
		csvX2 << "28";
	} else if (chrArr[i_goodpeak - 1] == "X") {
		csvX1 << "28.5";
		csvX2 << "29.5";
	} else if (chrArr[i_goodpeak - 1] == "Y") {
		csvX1 << "30";
		csvX2 << "31";
	} else { // different chromosomes written off grid boundaries
		csvX1 << "-5";
		csvX2 << "-5";
	}

	// band normalized from center to 50,000 bp XXX

	// ideoy1 position of band
	int i_mid = 0;
	for (int i = 0; i < i_goodpeak-1; i++) {
		i_mid = (chrEnd[i] - chrStart[i])/2;
		csvY1 << chrStart[i] + i_mid - 50000 << "\n";
	}
	i_mid = (chrEnd[i_goodpeak] - chrStart[i_goodpeak])/2;
	csvY1 << chrStart[i_goodpeak] + i_mid - 50000;

	// ideoy2 position of band
	i_mid = 0;
	for (int i = 0; i < i_goodpeak-1; i++) {
		i_mid = (chrEnd[i] - chrStart[i])/2;
		csvY2 << chrStart[i] + i_mid + 50000 << "\n";
	}
	i_mid = (chrEnd[i_goodpeak] - chrStart[i_goodpeak])/2;
	csvY2 << chrStart[i_goodpeak] + i_mid + 50000;



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

	rIdeo << "ch=data.frame(chx=c(0.5,2,3.5,5,6.5,8,9.5,11,12.5,14,15.5,17,18.5,20,21.5,23,24.5,26,27.5,29,30.5), t=c(\"1\",\"2\",\"3\",\"4\",\"5\",\"6\",\"7\",\"8\",\"9\",\"10\",\"11\",\"12\",\"13\",\"14\",\"15\",\"16\",\"17\",\"18\",\"19\",\"X\",\"Y\"))\n";
	rIdeo << "d=data.frame(x=c(0,1,1,0,0,1,1,0, 1.5,2.5,2.5,1.5,1.5,2.5,2.5,1.5, 3,4,4,3,3,4,4,3, 4.5,5.5,5.5,4.5,4.5,5.5,5.5,4.5, 6,7,7,6,6,7,7,6, 7.5,8.5,8.5,7.5,7.5,8.5,8.5,7.5, 9,10,10,9,9,10,10,9, 10.5,11.5,11.5,10.5,10.5,11.5,11.5,10.5, 12,13,13,12,12,13,13,12, 13.5,14.5,14.5,13.5,13.5,14.5,14.5,13.5, 15,16,16,15,15,16,16,15, 16.5,17.5,17.5,16.5,16.5,17.5,17.5,16.5, 18,19,19,18,18,19,19,18, 19.5,20.5,20.5,19.5,19.5,20.5,20.5,19.5, 21,22,22,21,21,22,22,21, 22.5,23.5,23.5,22.5,22.5,23.5,23.5,22.5, 24,25,25,24,24,25,25,24, 25.5,26.5,26.5,25.5,25.5,26.5,26.5,25.5, 27,28,28,27,27,28,28,27, 28.5,29.5,29.5,28.5,28.5,29.5,29.5,28.5, 30,31,31,30), y=c(110000,3000000,195471971,195471971,3000000,110000,0,0, 110000,3000000,182113224,182113224,3000000,110000,0,0, 110000,3000000,160039680,160039680,3000000,110000,0,0, 110000,3000000,156508116,156508116,3000000,110000,0,0, 110000,3000000,151834684,151834684,3000000,110000,0,0, 110000,3000000,149736546,149736546,3000000,110000,0,0, 110000,3000000,145441459,145441459,3000000,110000,0,0, 110000,3000000,129401213,129401213,3000000,110000,0,0, 110000,3000000,124595110,124595110,3000000,110000,0,0, 110000,3000000,130694993,130694993,3000000,110000,0,0, 110000,3000000,122082543,122082543,3000000,110000,0,0, 110000,3000000,120129022,120129022,3000000,110000,0,0, 110000,3000000,120421639,120421639,3000000,110000,0,0, 110000,3000000,124902244,124902244,3000000,110000,0,0, 110000,3000000,104043685,104043685,3000000,110000,0,0, 110000,3000000,98207768,98207768,3000000,110000,0,0, 110000,3000000,94987271,94987271,3000000,110000,0,0, 110000,3000000,90702639,90702639,3000000,110000,0,0, 110000,3000000,61431566,61431566,3000000,110000,0,0, 110000,3000000,171031299,171031299,3000000,110000,0,0, 0,0,91744698,91744698), t=c('1', '1', '1', '1', '1', '1', '1', '1', '2', '2', '2', '2', '2', '2', '2', '2', '3', '3', '3', '3', '3', '3', '3', '3', '4', '4', '4', '4', '4', '4', '4', '4', '5', '5', '5', '5', '5', '5', '5', '5', '6', '6', '6', '6', '6', '6', '6', '6', '7', '7', '7', '7', '7', '7', '7', '7', '8', '8', '8', '8', '8', '8', '8', '8', '9', '9', '9', '9', '9', '9', '9', '9', '10', '10', '10', '10', '10', '10', '10', '10', '11', '11', '11', '11', '11', '11', '11', '11', '12', '12', '12', '12', '12', '12', '12', '12', '13', '13', '13', '13', '13', '13', '13', '13', '14', '14', '14', '14', '14', '14', '14', '14', '15', '15', '15', '15', '15', '15', '15', '15', '16', '16', '16', '16', '16', '16', '16', '16', '17', '17', '17', '17', '17', '17', '17', '17', '18', '18', '18', '18', '18', '18', '18', '18', '19', '19', '19', '19', '19', '19', '19', '19', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'Y', 'Y', 'Y', 'Y'))\n";

	// this is same for all except the scale_x_continuous(limits = c(0, 32)) - 32 is variable

	rIdeo << "#y=c(cenS,cenE,chrE,chrE,cenE,cenS,0,0,\n";
	if (s_direction == "forward") { rIdeo << "ideForward <- ggplot() +\n"; }			
	if (s_direction == "reverse") { rIdeo << "ideReverse <- ggplot() +\n"; }			
	rIdeo << "  scale_x_continuous(limits = c(0, 32)) +\n";
	rIdeo << "  scale_y_continuous(limits = c(-5000000, 200000000)) +\n";
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
