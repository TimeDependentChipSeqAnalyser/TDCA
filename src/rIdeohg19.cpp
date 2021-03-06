#include <iostream>
#include <fstream>
#include <string>
//#include <regex> 
#include <vector>
#include <algorithm>
using namespace std;

std::string exec(const char* cmd);

// creates ideogram with trunover heatmap
void rIdeohg19(std::string s_rPltsName, int i_bamFiles, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, std::string s_direction)
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
			std::string s = dataArray[i][2*i_bamFiles + 5];		
			if ( (s != "NaN") && (s != "-") ) { i_goodpeak++; }	
		}
		if (s_direction == "reverse") {					
			std::string s = dataArray[i][2*i_bamFiles + 12];		
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
			std::string s_turnover = dataArray[i][2*i_bamFiles + 5];		
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
			std::string s_turnover = dataArray[i][2*i_bamFiles + 12];		
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


	rIdeo << "#hg19\n";
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
		} else if (chrArr[i] == "20") {
			csvX1 << "28.5\n";
			csvX2 << "29.5\n";
		} else if (chrArr[i] == "21") {
			csvX1 << "30\n";
			csvX2 << "31\n";
		} else if (chrArr[i] == "22") {
			csvX1 << "31.5\n";
			csvX2 << "32.5\n";
		} else if (chrArr[i] == "X") {
			csvX1 << "33\n";
			csvX2 << "34\n";
		} else if (chrArr[i] == "Y") {
			csvX1 << "34.5\n";
			csvX2 << "35.5\n";
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
	} else if (chrArr[i_goodpeak - 1] == "20") {
		csvX1 << "28.5";
		csvX2 << "29.5";
	} else if (chrArr[i_goodpeak - 1] == "21") {
		csvX1 << "30";
		csvX2 << "31";
	} else if (chrArr[i_goodpeak - 1] == "22") {
		csvX1 << "31.5";
		csvX2 << "32.5";
	} else if (chrArr[i_goodpeak - 1] == "X") {
		csvX1 << "33";
		csvX2 << "34\n";
	} else if (chrArr[i_goodpeak - 1] == "Y") {
		csvX1 << "34.5";
		csvX2 << "35.5";
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
			std::string s = dataArray[i][2*i_bamFiles + 5];					
			if ( (s != "NaN") && (s != "-") && (i_iterator == i_goodpeak-1) ) {
				csvInf << s;							
			}												
			if ( (s != "NaN") && (s != "-") && (i_iterator != i_goodpeak-1) ) {	
				csvInf << s << "\n";							
				i_iterator++;							
			}											
		}											
		if (s_direction == "reverse") {						
			std::string s = dataArray[i][2*i_bamFiles + 12];			
			if ( (s != "NaN") && (s != "-") && (i_iterator == i_goodpeak-1) ) {
				csvInf << s;						
			}	
			if ( (s != "NaN") && (s != "-") && (i_iterator != i_goodpeak-1) ) {
				csvInf << s << "\n";							
				i_iterator++;								
			}										
		}
	}


	// Draw all chromosomes XXX DOUBLE CHECK CENTROMERE POSITIONS 01/12/2017
	
	rIdeo << "ch=data.frame(chx=c(0.5,2,3.5,5,6.5,8,9.5,11,12.5,14,15.5,17,18.5,20,21.5,23,24.5,26,27.5,29,30.5,32,33.5,35), t=c(\"1\",\"2\",\"3\",\"4\",\"5\",\"6\",\"7\",\"8\",\"9\",\"10\",\"11\",\"12\",\"13\",\"14\",\"15\",\"16\",\"17\",\"18\",\"19\",\"20\",\"21\",\"22\",\"X\",\"Y\"))\n";
	rIdeo << "d=data.frame(x=c(0,1,1,0,0,1,1,0, 1.5,2.5,2.5,1.5,1.5,2.5,2.5,1.5, 3,4,4,3,3,4,4,3, 4.5,5.5,5.5,4.5,4.5,5.5,5.5,4.5, 6,7,7,6,6,7,7,6, 7.5,8.5,8.5,7.5,7.5,8.5,8.5,7.5, 9,10,10,9,9,10,10,9, 10.5,11.5,11.5,10.5,10.5,11.5,11.5,10.5, 12,13,13,12,12,13,13,12, 13.5,14.5,14.5,13.5,13.5,14.5,14.5,13.5, 15,16,16,15,15,16,16,15, 16.5,17.5,17.5,16.5,16.5,17.5,17.5,16.5, 18,19,19,18,18,19,19,18, 19.5,20.5,20.5,19.5,19.5,20.5,20.5,19.5, 21,22,22,21,21,22,22,21, 22.5,23.5,23.5,22.5,22.5,23.5,23.5,22.5, 24,25,25,24,24,25,25,24, 25.5,26.5,26.5,25.5,25.5,26.5,26.5,25.5, 27,28,28,27,27,28,28,27, 28.5,29.5,29.5,28.5,28.5,29.5,29.5,28.5, 30,31,31,30,30,31,31,30, 31.5,32.5,32.5,31.5,31.5,32.5,32.5,31.5, 33,34,34,33,33,34,34,33, 34.5,35.5,35.5,34.5,34.5,35.5,35.5,34.5), y=c(121535434,124535434,249250621,249250621,124535434,121535434,0,0, 92326171,95326171,243199373,243199373,95326171,92326171,0,0, 90504854,93504854,198022430,198022430,93504854,90504854,0,0, 49660117,52660117,191154276,191154276,52660117,49660117,0,0, 46405641,49405641,180915260,180915260,49405641,46405641,0,0, 58830166,61830166,171115067,171115067,61830166,58830166,0,0, 58054331,61054331,159138663,159138663,61054331,58054331,0,0, 43838887,46838887,146364022,146364022,46838887,43838887,0,0, 47367679,50367679,141213431,141213431,50367679,47367679,0,0, 39254935,42254935,135534747,135534747,42254935,39254935,0,0, 51644205,54644205,135006516,135006516,54644205,51644205,0,0, 34856694,37856694,133851895,133851895,37856694,34856694,0,0, 16000000,19000000,115169878,115169878,19000000,16000000,0,0, 16000000,19000000,107349540,107349540,19000000,16000000,0,0, 17000000,20000000,102531392,102531392,20000000,17000000,0,0, 35335801,38335801,90354753,90354753,38335801,35335801,0,0, 22263006,25263006,81195210,81195210,25263006,22263006,0,0, 15460898,18460898,78077248,78077248,18460898,15460898,0,0, 24681782,27681782,59128983,59128983,27681782,24681782,0,0, 26369569,29369569,63025520,63025520,29369569,26369569,0,0, 11288129,14288129,48129895,48129895,14288129,11288129,0,0, 13000000,16000000,51304566,51304566,16000000,13000000,0,0, 58632012,61632012,155270560,155270560,61632012,58632012,0,0, 10104553,13104553,59373566,59373566,13104553,10104553,0,0), t=c('1', '1', '1', '1', '1', '1', '1', '1', '2', '2', '2', '2', '2', '2', '2', '2', '3', '3', '3', '3', '3', '3', '3', '3', '4', '4', '4', '4', '4', '4', '4', '4', '5', '5', '5', '5', '5', '5', '5', '5', '6', '6', '6', '6', '6', '6', '6', '6', '7', '7', '7', '7', '7', '7', '7', '7', '8', '8', '8', '8', '8', '8', '8', '8', '9', '9', '9', '9', '9', '9', '9', '9', '10', '10', '10', '10', '10', '10', '10', '10', '11', '11', '11', '11', '11', '11', '11', '11', '12', '12', '12', '12', '12', '12', '12', '12', '13', '13', '13', '13', '13', '13', '13', '13', '14', '14', '14', '14', '14', '14', '14', '14', '15', '15', '15', '15', '15', '15', '15', '15', '16', '16', '16', '16', '16', '16', '16', '16', '17', '17', '17', '17', '17', '17', '17', '17', '18', '18', '18', '18', '18', '18', '18', '18', '19', '19', '19', '19', '19', '19', '19', '19', '20', '20', '20', '20', '20', '20', '20', '20', '21', '21', '21', '21', '21', '21', '21', '21', '22', '22', '22', '22', '22', '22', '22', '22', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'X', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y', 'Y'))\n";
	rIdeo << "#y=c(cenS,cenE,chrE,chrE,cenE,cenS,0,0,\n";
	if (s_direction == "forward") { rIdeo << "ideForward <- ggplot() +\n"; }			
	if (s_direction == "reverse") { rIdeo << "ideReverse <- ggplot() +\n"; }			
	rIdeo << "  scale_x_continuous(limits = c(0, 36)) +\n";
	rIdeo << "  scale_y_continuous(limits = c(-5000000, 250000000)) +\n";
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


















