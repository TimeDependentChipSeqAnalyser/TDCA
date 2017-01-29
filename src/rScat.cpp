#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>    // std::sort

using namespace std;

void rScat(std::string s_rPltsName, int i_bamFiles, int turnoverTimes[], int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, int i_undefRise, int i_undefFall, int i_valley, int i_hill, int i_fall, int i_rise, bool b_model, std::string s_direction)
{
	/** RECALL:
	dataArray[0][2*i_bamFiles + 3]  = "s_type_vec[i]"; 
	dataArray[0][2*i_bamFiles + 4]  = "s_drcFI_vec[i]"; 
	dataArray[0][2*i_bamFiles + 5]  = "s_drcFH_vec[i]";
	dataArray[0][2*i_bamFiles + 6]  = "s_drcFU_vec[i]"; 
	dataArray[0][2*i_bamFiles + 7]  = "s_drcFL_vec[i]";	
	dataArray[0][2*i_bamFiles + 8]  = "s_drcRI_vec[i]"; 
	dataArray[0][2*i_bamFiles + 9]  = "s_drcRH_vec[i]"; 
	dataArray[0][2*i_bamFiles + 10] = "s_drcRU_vec[i]"; 
	dataArray[0][2*i_bamFiles + 11] = "s_drcRL_vec[i]";
	**/

	// write R script
	ofstream rScat;
	rScat.open (s_rPltsName, std::ios::app);

	// write csv files 
	ofstream csvScatInf;
	ofstream csvScatHC;
	ofstream csvScatID;
	std::string s_csvScatInfName;
	std::string s_csvScatHCName;
	std::string s_csvScatIDName;

	if (s_direction == "forward") { s_csvScatInfName = "./" + s_rPltsName + "_data/fScatInf.csv"; }
	if (s_direction == "reverse") { s_csvScatInfName = "./" + s_rPltsName + "_data/rScatInf.csv"; }
	if (s_direction == "forward") { s_csvScatHCName = "./" + s_rPltsName + "_data/fScatHC.csv"; }
	if (s_direction == "reverse") { s_csvScatHCName = "./" + s_rPltsName + "_data/rScatHC.csv"; }
	if (s_direction == "forward") { s_csvScatIDName = "./" + s_rPltsName + "_data/fScatID.csv"; }
	if (s_direction == "reverse") { s_csvScatIDName = "./" + s_rPltsName + "_data/rScatID.csv"; }

	csvScatInf.open(s_csvScatInfName);
	csvScatHC.open(s_csvScatHCName);
	csvScatID.open(s_csvScatIDName);


	int i_forward = i_undefRise + i_rise + i_hill + i_valley;
	int i_reverse = i_undefFall + i_fall + i_valley + i_hill;

	int i_iterator = 0;
	// these are for point off the chart so that the legend remains accurate
	csvScatInf << "-1000\n-1000\n-1000\n-1000\n"; 

	for (int i = 1; i < i_peakNumber+1; i++) {	
		if (s_direction == "forward") {
			std::string s_i = dataArray[i][2*i_bamFiles + 5];
			if ( (s_i != "-") &&  (s_i != "NaN") ) {	 	// any rise
				if (i_iterator == i_forward-1) { 		// last data value
					csvScatInf << dataArray[i][2*i_bamFiles + 4];	
					i_iterator++;
				}
				if (i_iterator < i_forward-1) { 		// not last data value
					csvScatInf << dataArray[i][2*i_bamFiles + 4] << "\n";
					i_iterator++;
				}			
			}
		}
		if (s_direction == "reverse") {
			std::string s_i = dataArray[i][2*i_bamFiles + 9];
			if ( (s_i != "-") &&  (s_i != "NaN") ) {	 	// any fall
				if (i_iterator == i_reverse-1) { 		// last data value
					csvScatInf << dataArray[i][2*i_bamFiles + 8];
					i_iterator++;
				}
				if (i_iterator < i_reverse-1) { 		// not last data value
					csvScatInf << dataArray[i][2*i_bamFiles + 8] << "\n";
					i_iterator++;
				}			
			}
		}
	}

	i_iterator = 0;
	// these are for point off the chart so that the legend remains accurate
	csvScatHC << "-1000\n-1000\n-1000\n-1000\n";

	for (int i = 1; i < i_peakNumber+1; i++) {	
		if (s_direction == "forward") {
			std::string s_h = dataArray[i][2*i_bamFiles + 5];
			if ( (s_h != "-") &&  (s_h != "NaN") ) {	 	// any rise
				if (i_iterator == i_forward-1) { 		// last data value
					csvScatHC << dataArray[i][2*i_bamFiles + 5];
					i_iterator++;
				}
				if (i_iterator < i_forward-1) { 		// not last data value
					csvScatHC << dataArray[i][2*i_bamFiles + 5] << "\n";
					i_iterator++;
				}			
			}
		}
		if (s_direction == "reverse") {
			std::string s_h = dataArray[i][2*i_bamFiles + 9];
			if ( (s_h != "-") &&  (s_h != "NaN") ) {	 	// any fall
				if (i_iterator == i_reverse-1) { 		// last data value
					csvScatHC << dataArray[i][2*i_bamFiles + 9];
					i_iterator++;
				}
				if (i_iterator < i_reverse-1) { 		// not last data value
					csvScatHC << dataArray[i][2*i_bamFiles + 9] << "\n";
					i_iterator++;
				}			
			}
		}
	}
	std::vector<std::string> s_inflectionBinsTypeVec; 	// holds type (rise, fall, hill, etc.) of inflection data.
	for (int i = 1; i < i_peakNumber+1; i++) {
		if (s_direction == "forward") {
			std::string s = dataArray[i][2*i_bamFiles + 5];				
			if (dataArray[i][2*i_bamFiles + 3] == "rise") {	
				s_inflectionBinsTypeVec.push_back("rise"); }		
			if (dataArray[i][2*i_bamFiles + 3] == "hill") {					
				s_inflectionBinsTypeVec.push_back("hill"); }			
			if (dataArray[i][2*i_bamFiles + 3] == "valley") {				
				s_inflectionBinsTypeVec.push_back("valley"); }
			if ( (dataArray[i][2*i_bamFiles + 3] == "undefined") && (s != "-") &&  (s != "NaN") ) { 	
				s_inflectionBinsTypeVec.push_back("undefined rise"); }
		}
		if (s_direction == "reverse") {	
			std::string s = dataArray[i][2*i_bamFiles + 9];			
			if (dataArray[i][2*i_bamFiles + 3] == "fall") {	
				s_inflectionBinsTypeVec.push_back("fall"); }		
			if (dataArray[i][2*i_bamFiles + 3] == "hill") {					
				s_inflectionBinsTypeVec.push_back("hill"); }			
			if (dataArray[i][2*i_bamFiles + 3] == "valley") {				
				s_inflectionBinsTypeVec.push_back("valley"); }
			if ( (dataArray[i][2*i_bamFiles + 3] == "undefined") && (s != "-") &&  (s != "NaN") ) { 	
				s_inflectionBinsTypeVec.push_back("undefined fall"); }
		}
	}

	int i_types = s_inflectionBinsTypeVec.size();
	if (s_direction == "forward") { csvScatID << "undefined rise\nvalley\nhill\nrise\n"; }
	if (s_direction == "reverse") { csvScatID << "undefined fall\nvalley\nhill\nfall\n"; }

	for (int i = 0; i < i_types-1; i++)
		csvScatID << s_inflectionBinsTypeVec[i] << "\n";
	csvScatID << s_inflectionBinsTypeVec[i_types-1];

	if (s_direction == "forward") {
		rScat << "fScatInf <- read.csv(file=\"" << s_csvScatInfName << "\", header=FALSE)\n";
		rScat << "fScatHC <- read.csv(file=\"" << s_csvScatHCName << "\", header=FALSE)\n";
		rScat << "fScatID <- read.csv(file=\"" << s_csvScatIDName << "\", header=FALSE, stringsAsFactors=FALSE)\n";
		rScat << "fScatFrame<- data.frame(fScatHC,fScatInf,fScatID)\n"; 
	}
	if (s_direction == "reverse") {
		rScat << "rScatInf <- read.csv(file=\"" << s_csvScatInfName << "\", header=FALSE)\n";
		rScat << "rScatHC <- read.csv(file=\"" << s_csvScatHCName << "\", header=FALSE)\n";
		rScat << "rScatID <- read.csv(file=\"" << s_csvScatIDName << "\", header=FALSE, stringsAsFactors=FALSE)\n";
		rScat << "rScatFrame<- data.frame(rScatHC,rScatInf,rScatID)\n";
	}
	if (s_direction == "forward") {
		if (b_model) {
			rScat << "forwardScat <- ggplot(fScatFrame, aes(fScatFrame[,2],fScatFrame[,1],colour=fScatFrame[,3])) + \n"; 
			rScat << "scale_colour_manual(values=c(\"red\", \"blue\", \"chartreuse4\", \"darkgoldenrod1\")) + \n"; 
		} else {
			rScat << "forwardScat <- ggplot(fScatFrame, aes(fScatFrame[,2],fScatFrame[,1])) + \n"; 
		}
		rScat << "	ggtitle(\"TTI and Hills Coefficient of Signal Increase Loci\") +\n"; 
	}
	if (s_direction == "reverse") {
		if (b_model) {
			rScat << "reverseScat <- ggplot(rScatFrame, aes(rScatFrame[,2],rScatFrame[,1],colour=rScatFrame[,3])) + \n"; 
			rScat << "scale_colour_manual(values=c(\"blue\", \"red\", \"chartreuse4\", \"darkgoldenrod1\")) + \n"; 
		} else {
			rScat << "reverseScat <- ggplot(rScatFrame, aes(rScatFrame[,2],rScatFrame[,1])) + \n"; 
		}
		rScat << "	ggtitle(\"TTI and Hills Coefficient of Signal Decrease Loci\") +\n"; 
	}

	s_inflectionBinsTypeVec.clear();

	std::vector<double> d_hcVec; // NEW XXX

	// XXX STACKED BAR CHART
	// XXX inflectionBins = 0 to largest timepoint

	double d_hillCoefMax = 0;
	double d_bin = double(turnoverTimes[i_bamFiles-1]) / 25;

	// 4 rows, 25 columns, initialized to 0. 
	// rows are: 0) rise/fall 1) hill 2) valley 3) undefined
	int i_stackBarArr[4][25]={0};

	for (int i = 1; i < i_peakNumber+1; i++) {	
		if (s_direction == "forward") {
			std::string s = dataArray[i][2*i_bamFiles + 5]; 	 // forward hills coefficient
			if ( (s != "-") &&  (s != "NaN") ) {	 		 // any rise
				double d = stod(dataArray[i][2*i_bamFiles + 4]); // forward inflection point
				std::string s_t = dataArray[i][2*i_bamFiles + 3];
			
				if ( (d > 0) && (d <= (d_bin*1)) && (s_t == "rise") )	 		   { i_stackBarArr[0][0]++; }
				if ( (d > 0) && (d <= (d_bin*1)) && (s_t == "hill") )	 		   { i_stackBarArr[1][0]++; }
				if ( (d > 0) && (d <= (d_bin*1)) && (s_t == "valley") )		   { i_stackBarArr[2][0]++; }
				if ( (d > 0) && (d <= (d_bin*1)) && (s_t == "undefined") )		   { i_stackBarArr[3][0]++; }

				if ( (d > (d_bin*1)) && (d <= (d_bin*2)) && (s_t == "rise") )	   { i_stackBarArr[0][1]++; }
				if ( (d > (d_bin*1)) && (d <= (d_bin*2)) && (s_t == "hill") )	   { i_stackBarArr[1][1]++; }
				if ( (d > (d_bin*1)) && (d <= (d_bin*2)) && (s_t == "valley") )	   { i_stackBarArr[2][1]++; }
				if ( (d > (d_bin*1)) && (d <= (d_bin*2)) && (s_t == "undefined") )   { i_stackBarArr[3][1]++; }

				if ( (d > (d_bin*2)) && (d <= (d_bin*3)) && (s_t == "rise") )	   { i_stackBarArr[0][2]++; }
				if ( (d > (d_bin*2)) && (d <= (d_bin*3)) && (s_t == "hill") )	   { i_stackBarArr[1][2]++; }
				if ( (d > (d_bin*2)) && (d <= (d_bin*3)) && (s_t == "valley") )	   { i_stackBarArr[2][2]++; }
				if ( (d > (d_bin*2)) && (d <= (d_bin*3)) && (s_t == "undefined") )   { i_stackBarArr[3][2]++; }

				if ( (d > (d_bin*3)) && (d <= (d_bin*4)) && (s_t == "rise") )	   { i_stackBarArr[0][3]++; }
				if ( (d > (d_bin*3)) && (d <= (d_bin*4)) && (s_t == "hill") )	   { i_stackBarArr[1][3]++; }
				if ( (d > (d_bin*3)) && (d <= (d_bin*4)) && (s_t == "valley") )	   { i_stackBarArr[2][3]++; }
				if ( (d > (d_bin*3)) && (d <= (d_bin*4)) && (s_t == "undefined") )   { i_stackBarArr[3][3]++; }

				if ( (d > (d_bin*4)) && (d <= (d_bin*5)) && (s_t == "rise") )	   { i_stackBarArr[0][4]++; }
				if ( (d > (d_bin*4)) && (d <= (d_bin*5)) && (s_t == "hill") )	   { i_stackBarArr[1][4]++; }
				if ( (d > (d_bin*4)) && (d <= (d_bin*5)) && (s_t == "valley") )	   { i_stackBarArr[2][4]++; }
				if ( (d > (d_bin*4)) && (d <= (d_bin*5)) && (s_t == "undefined") )   { i_stackBarArr[3][4]++; }

				if ( (d > (d_bin*5)) && (d <= (d_bin*6)) && (s_t == "rise") )	   { i_stackBarArr[0][5]++; }
				if ( (d > (d_bin*5)) && (d <= (d_bin*6)) && (s_t == "hill") )	   { i_stackBarArr[1][5]++; }
				if ( (d > (d_bin*5)) && (d <= (d_bin*6)) && (s_t == "valley") )	   { i_stackBarArr[2][5]++; }
				if ( (d > (d_bin*5)) && (d <= (d_bin*6)) && (s_t == "undefined") )   { i_stackBarArr[3][5]++; }

				if ( (d > (d_bin*6)) && (d <= (d_bin*7)) && (s_t == "rise") )	   { i_stackBarArr[0][6]++; }
				if ( (d > (d_bin*6)) && (d <= (d_bin*7)) && (s_t == "hill") )	   { i_stackBarArr[1][6]++; }
				if ( (d > (d_bin*6)) && (d <= (d_bin*7)) && (s_t == "valley") )	   { i_stackBarArr[2][6]++; }
				if ( (d > (d_bin*6)) && (d <= (d_bin*7)) && (s_t == "undefined") )   { i_stackBarArr[3][6]++; }

				if ( (d > (d_bin*7)) && (d <= (d_bin*8)) && (s_t == "rise") )	   { i_stackBarArr[0][7]++; }
				if ( (d > (d_bin*7)) && (d <= (d_bin*8)) && (s_t == "hill") )	   { i_stackBarArr[1][7]++; }
				if ( (d > (d_bin*7)) && (d <= (d_bin*8)) && (s_t == "valley") )	   { i_stackBarArr[2][7]++; }
				if ( (d > (d_bin*7)) && (d <= (d_bin*8)) && (s_t == "undefined") )   { i_stackBarArr[3][7]++; }

				if ( (d > (d_bin*8)) && (d <= (d_bin*9)) && (s_t == "rise") )	   { i_stackBarArr[0][8]++; }
				if ( (d > (d_bin*8)) && (d <= (d_bin*9)) && (s_t == "hill") )	   { i_stackBarArr[1][8]++; }
				if ( (d > (d_bin*8)) && (d <= (d_bin*9)) && (s_t == "valley") )	   { i_stackBarArr[2][8]++; }
				if ( (d > (d_bin*8)) && (d <= (d_bin*9)) && (s_t == "undefined") )   { i_stackBarArr[3][8]++; }

				if ( (d > (d_bin*9)) && (d <= (d_bin*10)) && (s_t == "rise") )	   { i_stackBarArr[0][9]++; }
				if ( (d > (d_bin*9)) && (d <= (d_bin*10)) && (s_t == "hill") )	   { i_stackBarArr[1][9]++; }
				if ( (d > (d_bin*9)) && (d <= (d_bin*10)) && (s_t == "valley") )	   { i_stackBarArr[2][9]++; }
				if ( (d > (d_bin*9)) && (d <= (d_bin*10)) && (s_t == "undefined") )  { i_stackBarArr[3][9]++; }

				if ( (d > (d_bin*10)) && (d <= (d_bin*11)) && (s_t == "rise") )	   { i_stackBarArr[0][10]++; }
				if ( (d > (d_bin*10)) && (d <= (d_bin*11)) && (s_t == "hill") )	   { i_stackBarArr[1][10]++; }
				if ( (d > (d_bin*10)) && (d <= (d_bin*11)) && (s_t == "valley") )	   { i_stackBarArr[2][10]++; }
				if ( (d > (d_bin*10)) && (d <= (d_bin*11)) && (s_t == "undefined") ) { i_stackBarArr[3][10]++; }

				if ( (d > (d_bin*11)) && (d <= (d_bin*12)) && (s_t == "rise") )	   { i_stackBarArr[0][11]++; }
				if ( (d > (d_bin*11)) && (d <= (d_bin*12)) && (s_t == "hill") )	   { i_stackBarArr[1][11]++; }
				if ( (d > (d_bin*11)) && (d <= (d_bin*12)) && (s_t == "valley") )	   { i_stackBarArr[2][11]++; }
				if ( (d > (d_bin*11)) && (d <= (d_bin*12)) && (s_t == "undefined") ) { i_stackBarArr[3][11]++; }

				if ( (d > (d_bin*12)) && (d <= (d_bin*13)) && (s_t == "rise") )	   { i_stackBarArr[0][12]++; }
				if ( (d > (d_bin*12)) && (d <= (d_bin*13)) && (s_t == "hill") )	   { i_stackBarArr[1][12]++; }
				if ( (d > (d_bin*12)) && (d <= (d_bin*13)) && (s_t == "valley") )	   { i_stackBarArr[2][12]++; }
				if ( (d > (d_bin*12)) && (d <= (d_bin*13)) && (s_t == "undefined") ) { i_stackBarArr[3][12]++; }

				if ( (d > (d_bin*13)) && (d <= (d_bin*14)) && (s_t == "rise") )	   { i_stackBarArr[0][13]++; }
				if ( (d > (d_bin*13)) && (d <= (d_bin*14)) && (s_t == "hill") )	   { i_stackBarArr[1][13]++; }
				if ( (d > (d_bin*13)) && (d <= (d_bin*14)) && (s_t == "valley") )	   { i_stackBarArr[2][13]++; }
				if ( (d > (d_bin*13)) && (d <= (d_bin*14)) && (s_t == "undefined") ) { i_stackBarArr[3][13]++; }

				if ( (d > (d_bin*14)) && (d <= (d_bin*15)) && (s_t == "rise") )	   { i_stackBarArr[0][14]++; }
				if ( (d > (d_bin*14)) && (d <= (d_bin*15)) && (s_t == "hill") )	   { i_stackBarArr[1][14]++; }
				if ( (d > (d_bin*14)) && (d <= (d_bin*15)) && (s_t == "valley") )	   { i_stackBarArr[2][14]++; }
				if ( (d > (d_bin*14)) && (d <= (d_bin*15)) && (s_t == "undefined") ) { i_stackBarArr[3][14]++; }

				if ( (d > (d_bin*15)) && (d <= (d_bin*16)) && (s_t == "rise") )	   { i_stackBarArr[0][15]++; }
				if ( (d > (d_bin*15)) && (d <= (d_bin*16)) && (s_t == "hill") )	   { i_stackBarArr[1][15]++; }
				if ( (d > (d_bin*15)) && (d <= (d_bin*16)) && (s_t == "valley") )	   { i_stackBarArr[2][15]++; }
				if ( (d > (d_bin*15)) && (d <= (d_bin*16)) && (s_t == "undefined") ) { i_stackBarArr[3][15]++; }

				if ( (d > (d_bin*16)) && (d <= (d_bin*17)) && (s_t == "rise") )	   { i_stackBarArr[0][16]++; }
				if ( (d > (d_bin*16)) && (d <= (d_bin*17)) && (s_t == "hill") )	   { i_stackBarArr[1][16]++; }
				if ( (d > (d_bin*16)) && (d <= (d_bin*17)) && (s_t == "valley") )	   { i_stackBarArr[2][16]++; }
				if ( (d > (d_bin*16)) && (d <= (d_bin*17)) && (s_t == "undefined") ) { i_stackBarArr[3][16]++; }

				if ( (d > (d_bin*17)) && (d <= (d_bin*18)) && (s_t == "rise") )	   { i_stackBarArr[0][17]++; }
				if ( (d > (d_bin*17)) && (d <= (d_bin*18)) && (s_t == "hill") )	   { i_stackBarArr[1][17]++; }
				if ( (d > (d_bin*17)) && (d <= (d_bin*18)) && (s_t == "valley") )	   { i_stackBarArr[2][17]++; }
				if ( (d > (d_bin*17)) && (d <= (d_bin*18)) && (s_t == "undefined") ) { i_stackBarArr[3][17]++; }

				if ( (d > (d_bin*18)) && (d <= (d_bin*19)) && (s_t == "rise") )	   { i_stackBarArr[0][18]++; }
				if ( (d > (d_bin*18)) && (d <= (d_bin*19)) && (s_t == "hill") )	   { i_stackBarArr[1][18]++; }
				if ( (d > (d_bin*18)) && (d <= (d_bin*19)) && (s_t == "valley") )	   { i_stackBarArr[2][18]++; }
				if ( (d > (d_bin*18)) && (d <= (d_bin*19)) && (s_t == "undefined") ) { i_stackBarArr[3][18]++; }

				if ( (d > (d_bin*19)) && (d <= (d_bin*20)) && (s_t == "rise") )	   { i_stackBarArr[0][19]++; }
				if ( (d > (d_bin*19)) && (d <= (d_bin*20)) && (s_t == "hill") )	   { i_stackBarArr[1][19]++; }
				if ( (d > (d_bin*19)) && (d <= (d_bin*20)) && (s_t == "valley") )	   { i_stackBarArr[2][19]++; }
				if ( (d > (d_bin*19)) && (d <= (d_bin*20)) && (s_t == "undefined") ) { i_stackBarArr[3][19]++; }

				if ( (d > (d_bin*20)) && (d <= (d_bin*21)) && (s_t == "rise") )	   { i_stackBarArr[0][20]++; }
				if ( (d > (d_bin*20)) && (d <= (d_bin*21)) && (s_t == "hill") )	   { i_stackBarArr[1][20]++; }
				if ( (d > (d_bin*20)) && (d <= (d_bin*21)) && (s_t == "valley") )	   { i_stackBarArr[2][20]++; }
				if ( (d > (d_bin*20)) && (d <= (d_bin*21)) && (s_t == "undefined") ) { i_stackBarArr[3][20]++; }

				if ( (d > (d_bin*21)) && (d <= (d_bin*22)) && (s_t == "rise") )	   { i_stackBarArr[0][21]++; }
				if ( (d > (d_bin*21)) && (d <= (d_bin*22)) && (s_t == "hill") )	   { i_stackBarArr[1][21]++; }
				if ( (d > (d_bin*21)) && (d <= (d_bin*22)) && (s_t == "valley") )	   { i_stackBarArr[2][21]++; }
				if ( (d > (d_bin*21)) && (d <= (d_bin*22)) && (s_t == "undefined") ) { i_stackBarArr[3][21]++; }

				if ( (d > (d_bin*22)) && (d <= (d_bin*23)) && (s_t == "rise") )	   { i_stackBarArr[0][22]++; }
				if ( (d > (d_bin*22)) && (d <= (d_bin*23)) && (s_t == "hill") )	   { i_stackBarArr[1][22]++; }
				if ( (d > (d_bin*22)) && (d <= (d_bin*23)) && (s_t == "valley") )	   { i_stackBarArr[2][22]++; }
				if ( (d > (d_bin*22)) && (d <= (d_bin*23)) && (s_t == "undefined") ) { i_stackBarArr[3][22]++; }

				if ( (d > (d_bin*23)) && (d <= (d_bin*24)) && (s_t == "rise") )	   { i_stackBarArr[0][23]++; }
				if ( (d > (d_bin*23)) && (d <= (d_bin*24)) && (s_t == "hill") )	   { i_stackBarArr[1][23]++; }
				if ( (d > (d_bin*23)) && (d <= (d_bin*24)) && (s_t == "valley") )	   { i_stackBarArr[2][23]++; }
				if ( (d > (d_bin*23)) && (d <= (d_bin*24)) && (s_t == "undefined") ) { i_stackBarArr[3][23]++; }

				if ( (d > (d_bin*24)) && (d <= (d_bin*25)) && (s_t == "rise") )	   { i_stackBarArr[0][24]++; }
				if ( (d > (d_bin*24)) && (d <= (d_bin*25)) && (s_t == "hill") )	   { i_stackBarArr[1][24]++; }
				if ( (d > (d_bin*24)) && (d <= (d_bin*25)) && (s_t == "valley") )	   { i_stackBarArr[2][24]++; }
				if ( (d > (d_bin*24)) && (d <= (d_bin*25)) && (s_t == "undefined") ) { i_stackBarArr[3][24]++; }
			
				// simultaneously check hills coefficient
				if ( (d > 0) && (d <= (d_bin*25)) ) {
					double d_co = stod(dataArray[i][2*i_bamFiles + 5]);
					if (d_co < d_hillCoefMax) { d_hillCoefMax = d_co; }
					d_hcVec.push_back(d_co); // NEW XXX				
				}
			
			}
		}

		if (s_direction == "reverse") {
			std::string s = dataArray[i][2*i_bamFiles + 9]; 	 // reverse hills coefficien
			if ( (s != "-") &&  (s != "NaN") ) {	 		 // any fall
				double d = stod(dataArray[i][2*i_bamFiles + 8]); // reverse inflection point
				std::string s_t = dataArray[i][2*i_bamFiles + 3];
			
				if ( (d > 0) && (d <= (d_bin*1)) && (s_t == "fall") )	 		   { i_stackBarArr[0][0]++; }
				if ( (d > 0) && (d <= (d_bin*1)) && (s_t == "hill") )	 		   { i_stackBarArr[1][0]++; }
				if ( (d > 0) && (d <= (d_bin*1)) && (s_t == "valley") )		   { i_stackBarArr[2][0]++; }
				if ( (d > 0) && (d <= (d_bin*1)) && (s_t == "undefined") )		   { i_stackBarArr[3][0]++; }

				if ( (d > (d_bin*1)) && (d <= (d_bin*2)) && (s_t == "fall") )	   { i_stackBarArr[0][1]++; }
				if ( (d > (d_bin*1)) && (d <= (d_bin*2)) && (s_t == "hill") )	   { i_stackBarArr[1][1]++; }
				if ( (d > (d_bin*1)) && (d <= (d_bin*2)) && (s_t == "valley") )	   { i_stackBarArr[2][1]++; }
				if ( (d > (d_bin*1)) && (d <= (d_bin*2)) && (s_t == "undefined") )   { i_stackBarArr[3][1]++; }

				if ( (d > (d_bin*2)) && (d <= (d_bin*3)) && (s_t == "fall") )	   { i_stackBarArr[0][2]++; }
				if ( (d > (d_bin*2)) && (d <= (d_bin*3)) && (s_t == "hill") )	   { i_stackBarArr[1][2]++; }
				if ( (d > (d_bin*2)) && (d <= (d_bin*3)) && (s_t == "valley") )	   { i_stackBarArr[2][2]++; }
				if ( (d > (d_bin*2)) && (d <= (d_bin*3)) && (s_t == "undefined") )   { i_stackBarArr[3][2]++; }

				if ( (d > (d_bin*3)) && (d <= (d_bin*4)) && (s_t == "fall") )	   { i_stackBarArr[0][3]++; }
				if ( (d > (d_bin*3)) && (d <= (d_bin*4)) && (s_t == "hill") )	   { i_stackBarArr[1][3]++; }
				if ( (d > (d_bin*3)) && (d <= (d_bin*4)) && (s_t == "valley") )	   { i_stackBarArr[2][3]++; }
				if ( (d > (d_bin*3)) && (d <= (d_bin*4)) && (s_t == "undefined") )   { i_stackBarArr[3][3]++; }

				if ( (d > (d_bin*4)) && (d <= (d_bin*5)) && (s_t == "fall") )	   { i_stackBarArr[0][4]++; }
				if ( (d > (d_bin*4)) && (d <= (d_bin*5)) && (s_t == "hill") )	   { i_stackBarArr[1][4]++; }
				if ( (d > (d_bin*4)) && (d <= (d_bin*5)) && (s_t == "valley") )	   { i_stackBarArr[2][4]++; }
				if ( (d > (d_bin*4)) && (d <= (d_bin*5)) && (s_t == "undefined") )   { i_stackBarArr[3][4]++; }

				if ( (d > (d_bin*5)) && (d <= (d_bin*6)) && (s_t == "fall") )	   { i_stackBarArr[0][5]++; }
				if ( (d > (d_bin*5)) && (d <= (d_bin*6)) && (s_t == "hill") )	   { i_stackBarArr[1][5]++; }
				if ( (d > (d_bin*5)) && (d <= (d_bin*6)) && (s_t == "valley") )	   { i_stackBarArr[2][5]++; }
				if ( (d > (d_bin*5)) && (d <= (d_bin*6)) && (s_t == "undefined") )   { i_stackBarArr[3][5]++; }

				if ( (d > (d_bin*6)) && (d <= (d_bin*7)) && (s_t == "fall") )	   { i_stackBarArr[0][6]++; }
				if ( (d > (d_bin*6)) && (d <= (d_bin*7)) && (s_t == "hill") )	   { i_stackBarArr[1][6]++; }
				if ( (d > (d_bin*6)) && (d <= (d_bin*7)) && (s_t == "valley") )	   { i_stackBarArr[2][6]++; }
				if ( (d > (d_bin*6)) && (d <= (d_bin*7)) && (s_t == "undefined") )   { i_stackBarArr[3][6]++; }

				if ( (d > (d_bin*7)) && (d <= (d_bin*8)) && (s_t == "fall") )	   { i_stackBarArr[0][7]++; }
				if ( (d > (d_bin*7)) && (d <= (d_bin*8)) && (s_t == "hill") )	   { i_stackBarArr[1][7]++; }
				if ( (d > (d_bin*7)) && (d <= (d_bin*8)) && (s_t == "valley") )	   { i_stackBarArr[2][7]++; }
				if ( (d > (d_bin*7)) && (d <= (d_bin*8)) && (s_t == "undefined") )   { i_stackBarArr[3][7]++; }

				if ( (d > (d_bin*8)) && (d <= (d_bin*9)) && (s_t == "fall") )	   { i_stackBarArr[0][8]++; }
				if ( (d > (d_bin*8)) && (d <= (d_bin*9)) && (s_t == "hill") )	   { i_stackBarArr[1][8]++; }
				if ( (d > (d_bin*8)) && (d <= (d_bin*9)) && (s_t == "valley") )	   { i_stackBarArr[2][8]++; }
				if ( (d > (d_bin*8)) && (d <= (d_bin*9)) && (s_t == "undefined") )   { i_stackBarArr[3][8]++; }

				if ( (d > (d_bin*9)) && (d <= (d_bin*10)) && (s_t == "fall") )	   { i_stackBarArr[0][9]++; }
				if ( (d > (d_bin*9)) && (d <= (d_bin*10)) && (s_t == "hill") )	   { i_stackBarArr[1][9]++; }
				if ( (d > (d_bin*9)) && (d <= (d_bin*10)) && (s_t == "valley") )	   { i_stackBarArr[2][9]++; }
				if ( (d > (d_bin*9)) && (d <= (d_bin*10)) && (s_t == "undefined") )  { i_stackBarArr[3][9]++; }

				if ( (d > (d_bin*10)) && (d <= (d_bin*11)) && (s_t == "fall") )	   { i_stackBarArr[0][10]++; }
				if ( (d > (d_bin*10)) && (d <= (d_bin*11)) && (s_t == "hill") )	   { i_stackBarArr[1][10]++; }
				if ( (d > (d_bin*10)) && (d <= (d_bin*11)) && (s_t == "valley") )	   { i_stackBarArr[2][10]++; }
				if ( (d > (d_bin*10)) && (d <= (d_bin*11)) && (s_t == "undefined") ) { i_stackBarArr[3][10]++; }

				if ( (d > (d_bin*11)) && (d <= (d_bin*12)) && (s_t == "fall") )	   { i_stackBarArr[0][11]++; }
				if ( (d > (d_bin*11)) && (d <= (d_bin*12)) && (s_t == "hill") )	   { i_stackBarArr[1][11]++; }
				if ( (d > (d_bin*11)) && (d <= (d_bin*12)) && (s_t == "valley") )	   { i_stackBarArr[2][11]++; }
				if ( (d > (d_bin*11)) && (d <= (d_bin*12)) && (s_t == "undefined") ) { i_stackBarArr[3][11]++; }

				if ( (d > (d_bin*12)) && (d <= (d_bin*13)) && (s_t == "fall") )	   { i_stackBarArr[0][12]++; }
				if ( (d > (d_bin*12)) && (d <= (d_bin*13)) && (s_t == "hill") )	   { i_stackBarArr[1][12]++; }
				if ( (d > (d_bin*12)) && (d <= (d_bin*13)) && (s_t == "valley") )	   { i_stackBarArr[2][12]++; }
				if ( (d > (d_bin*12)) && (d <= (d_bin*13)) && (s_t == "undefined") ) { i_stackBarArr[3][12]++; }

				if ( (d > (d_bin*13)) && (d <= (d_bin*14)) && (s_t == "fall") )	   { i_stackBarArr[0][13]++; }
				if ( (d > (d_bin*13)) && (d <= (d_bin*14)) && (s_t == "hill") )	   { i_stackBarArr[1][13]++; }
				if ( (d > (d_bin*13)) && (d <= (d_bin*14)) && (s_t == "valley") )	   { i_stackBarArr[2][13]++; }
				if ( (d > (d_bin*13)) && (d <= (d_bin*14)) && (s_t == "undefined") ) { i_stackBarArr[3][13]++; }

				if ( (d > (d_bin*14)) && (d <= (d_bin*15)) && (s_t == "fall") )	   { i_stackBarArr[0][14]++; }
				if ( (d > (d_bin*14)) && (d <= (d_bin*15)) && (s_t == "hill") )	   { i_stackBarArr[1][14]++; }
				if ( (d > (d_bin*14)) && (d <= (d_bin*15)) && (s_t == "valley") )	   { i_stackBarArr[2][14]++; }
				if ( (d > (d_bin*14)) && (d <= (d_bin*15)) && (s_t == "undefined") ) { i_stackBarArr[3][14]++; }

				if ( (d > (d_bin*15)) && (d <= (d_bin*16)) && (s_t == "fall") )	   { i_stackBarArr[0][15]++; }
				if ( (d > (d_bin*15)) && (d <= (d_bin*16)) && (s_t == "hill") )	   { i_stackBarArr[1][15]++; }
				if ( (d > (d_bin*15)) && (d <= (d_bin*16)) && (s_t == "valley") )	   { i_stackBarArr[2][15]++; }
				if ( (d > (d_bin*15)) && (d <= (d_bin*16)) && (s_t == "undefined") ) { i_stackBarArr[3][15]++; }

				if ( (d > (d_bin*16)) && (d <= (d_bin*17)) && (s_t == "fall") )	   { i_stackBarArr[0][16]++; }
				if ( (d > (d_bin*16)) && (d <= (d_bin*17)) && (s_t == "hill") )	   { i_stackBarArr[1][16]++; }
				if ( (d > (d_bin*16)) && (d <= (d_bin*17)) && (s_t == "valley") )	   { i_stackBarArr[2][16]++; }
				if ( (d > (d_bin*16)) && (d <= (d_bin*17)) && (s_t == "undefined") ) { i_stackBarArr[3][16]++; }

				if ( (d > (d_bin*17)) && (d <= (d_bin*18)) && (s_t == "fall") )	   { i_stackBarArr[0][17]++; }
				if ( (d > (d_bin*17)) && (d <= (d_bin*18)) && (s_t == "hill") )	   { i_stackBarArr[1][17]++; }
				if ( (d > (d_bin*17)) && (d <= (d_bin*18)) && (s_t == "valley") )	   { i_stackBarArr[2][17]++; }
				if ( (d > (d_bin*17)) && (d <= (d_bin*18)) && (s_t == "undefined") ) { i_stackBarArr[3][17]++; }

				if ( (d > (d_bin*18)) && (d <= (d_bin*19)) && (s_t == "fall") )	   { i_stackBarArr[0][18]++; }
				if ( (d > (d_bin*18)) && (d <= (d_bin*19)) && (s_t == "hill") )	   { i_stackBarArr[1][18]++; }
				if ( (d > (d_bin*18)) && (d <= (d_bin*19)) && (s_t == "valley") )	   { i_stackBarArr[2][18]++; }
				if ( (d > (d_bin*18)) && (d <= (d_bin*19)) && (s_t == "undefined") ) { i_stackBarArr[3][18]++; }

				if ( (d > (d_bin*19)) && (d <= (d_bin*20)) && (s_t == "fall") )	   { i_stackBarArr[0][19]++; }
				if ( (d > (d_bin*19)) && (d <= (d_bin*20)) && (s_t == "hill") )	   { i_stackBarArr[1][19]++; }
				if ( (d > (d_bin*19)) && (d <= (d_bin*20)) && (s_t == "valley") )	   { i_stackBarArr[2][19]++; }
				if ( (d > (d_bin*19)) && (d <= (d_bin*20)) && (s_t == "undefined") ) { i_stackBarArr[3][19]++; }

				if ( (d > (d_bin*20)) && (d <= (d_bin*21)) && (s_t == "fall") )	   { i_stackBarArr[0][20]++; }
				if ( (d > (d_bin*20)) && (d <= (d_bin*21)) && (s_t == "hill") )	   { i_stackBarArr[1][20]++; }
				if ( (d > (d_bin*20)) && (d <= (d_bin*21)) && (s_t == "valley") )	   { i_stackBarArr[2][20]++; }
				if ( (d > (d_bin*20)) && (d <= (d_bin*21)) && (s_t == "undefined") ) { i_stackBarArr[3][20]++; }

				if ( (d > (d_bin*21)) && (d <= (d_bin*22)) && (s_t == "fall") )	   { i_stackBarArr[0][21]++; }
				if ( (d > (d_bin*21)) && (d <= (d_bin*22)) && (s_t == "hill") )	   { i_stackBarArr[1][21]++; }
				if ( (d > (d_bin*21)) && (d <= (d_bin*22)) && (s_t == "valley") )	   { i_stackBarArr[2][21]++; }
				if ( (d > (d_bin*21)) && (d <= (d_bin*22)) && (s_t == "undefined") ) { i_stackBarArr[3][21]++; }

				if ( (d > (d_bin*22)) && (d <= (d_bin*23)) && (s_t == "fall") )	   { i_stackBarArr[0][22]++; }
				if ( (d > (d_bin*22)) && (d <= (d_bin*23)) && (s_t == "hill") )	   { i_stackBarArr[1][22]++; }
				if ( (d > (d_bin*22)) && (d <= (d_bin*23)) && (s_t == "valley") )	   { i_stackBarArr[2][22]++; }
				if ( (d > (d_bin*22)) && (d <= (d_bin*23)) && (s_t == "undefined") ) { i_stackBarArr[3][22]++; }

				if ( (d > (d_bin*23)) && (d <= (d_bin*24)) && (s_t == "fall") )	   { i_stackBarArr[0][23]++; }
				if ( (d > (d_bin*23)) && (d <= (d_bin*24)) && (s_t == "hill") )	   { i_stackBarArr[1][23]++; }
				if ( (d > (d_bin*23)) && (d <= (d_bin*24)) && (s_t == "valley") )	   { i_stackBarArr[2][23]++; }
				if ( (d > (d_bin*23)) && (d <= (d_bin*24)) && (s_t == "undefined") ) { i_stackBarArr[3][23]++; }

				if ( (d > (d_bin*24)) && (d <= (d_bin*25)) && (s_t == "fall") )	   { i_stackBarArr[0][24]++; }
				if ( (d > (d_bin*24)) && (d <= (d_bin*25)) && (s_t == "hill") )	   { i_stackBarArr[1][24]++; }
				if ( (d > (d_bin*24)) && (d <= (d_bin*25)) && (s_t == "valley") )	   { i_stackBarArr[2][24]++; }
				if ( (d > (d_bin*24)) && (d <= (d_bin*25)) && (s_t == "undefined") ) { i_stackBarArr[3][24]++; }

				// simultaneously check hills coefficient
				if ( (d > 0) && (d <= (d_bin*25)) ) {
					double d_co = stod(dataArray[i][2*i_bamFiles + 9]);
					if (d_co > d_hillCoefMax) { d_hillCoefMax = d_co; }
					d_hcVec.push_back(d_co); // NEW XXX								
				}
			}
		}
	}

	// X axis limit to 90th percentile of data // NEW XXX
	std::sort (d_hcVec.begin(), d_hcVec.end()); // NEW XXX  
	int i_index90per = int(d_hcVec.size()*0.9); // NEW XXX
	int i_index10per = int(d_hcVec.size()*0.1); // NEW XXX

	// STACKED BAR OF HILL COEFICIENT 
	//double d_hcBin = d_hillCoefMax / 25;
	double d_hcBin; // NEW XXX
	if (s_direction == "forward") { 
		if (d_hcVec.size() < 10) { d_hcBin = -0.05; }	
		else { d_hcBin = d_hcVec[i_index10per] / 25; }
	}
	if (s_direction == "reverse") { 
		if (d_hcVec.size() < 10) { d_hcBin = 0.05; }
		else { d_hcBin = d_hcVec[i_index90per] / 25; }
	}

	rScat << "	geom_point(size = 0.3, alpha = 0.25, na.rm = T) +\n"; 
	rScat << "	guides(colour = guide_legend(override.aes = list(size = 1, alpha = 1))) +\n"; 
	rScat << "	labs(x = \"TTI (time)\",size = 10) +\n"; 
	rScat << "	labs(y = \"Hill's Coefficient\",size = 10) +\n"; 
	// max inflection point axis limit at 2*largest timepoint
	rScat << "	scale_x_continuous(expand = c(0,0), limits = c(0, " << turnoverTimes[i_bamFiles-1] << ")) +\n"; 
	if (s_direction == "forward") {
		if (d_hcVec.size() < 10) { rScat << "	scale_y_continuous(expand = c(0,0), limits = c(-1, 0)) +\n"; }
		else {rScat << "	scale_y_continuous(expand = c(0,0), limits = c(" << d_hcVec[i_index10per] << ", 0)) +\n"; }
	}
	if (s_direction == "reverse") {
		if (d_hcVec.size() < 10) { rScat << "	scale_y_continuous(expand = c(0,0), limits = c(0, 1)) +\n"; }
		else { rScat << "	scale_y_continuous(expand = c(0,0), limits = c(0, " << d_hcVec[i_index90per] << ")) +\n"; }
	}
	rScat << "	theme(\n"; 
	rScat << "		axis.text = element_text(size = 8, colour = \"black\"),\n"; 
	rScat << "		legend.key = element_rect(fill = \"white\"),\n"; 
	rScat << "		legend.background = element_rect(fill = \"white\"),\n"; 
	rScat << "		panel.grid.major = element_line(colour = \"white\"),\n"; 
	rScat << "		legend.position = \"\",\n"; 
	rScat << "		legend.title = element_text(size = 0),\n"; 
	rScat << "		panel.grid.minor = element_blank(),\n"; 
	rScat << "		panel.background = element_rect(fill = \"white\"),\n"; 
	rScat << "		axis.line.x = element_line(colour = \"black\"),\n"; 
	rScat << "		axis.line.y = element_line(colour = \"black\")\n"; 
	rScat << "	)\n"; 

	// rows are: 0) rise/fall 1) hill 2) valley 3) undefined
	rScat << "bin <- c(";
	for (int j = 0; j < 3; j++) {	
		for (int i = 0; i < 25; i++) {	
			rScat << (d_bin*i) + d_bin/2 << ", ";
		}
	}
	for (int i = 0; i < 24; i++) {	
		rScat << (d_bin*i) + d_bin/2 << ", ";
	}
	rScat << (d_bin*24) + d_bin/2 << ")\n";

	// print count
	rScat << "count <- c(";
	for (int i = 0; i < 3; i++) {	
		for (int j = 0; j < 25; j++) {	
			rScat << i_stackBarArr[i][j] << ", ";
		}
	}
	for (int i = 0; i < 24; i++) {	
		rScat << i_stackBarArr[3][i] << ", ";
	}
	rScat << i_stackBarArr[3][24] << ")\n";
	
	// print variable 
	rScat << "variable <- c(\n";
	if (s_direction == "forward") {
		for (int i = 0; i < 25; i++) {	
			rScat << "\"rise\", ";
		}
		for (int i = 0; i < 25; i++) {	
			rScat << "\"hill\", ";
		}
		for (int i = 0; i < 25; i++) {	
			rScat << "\"valley\", ";
		}
		for (int i = 0; i < 24; i++) {	
			rScat << "\"undefined rise\", ";
		}
		rScat << "\"undefined rise\")\n";
	}
	if (s_direction == "reverse") {
		for (int i = 0; i < 25; i++) {	
			rScat << "\"fall\", ";
		}
		for (int i = 0; i < 25; i++) {	
			rScat << "\"hill\", ";
		}
		for (int i = 0; i < 25; i++) {	
			rScat << "\"valley\", ";
		}
		for (int i = 0; i < 24; i++) {	
			rScat << "\"undefined fall\", ";
		}
		rScat << "\"undefined fall\")\n";
	}

	if (s_direction == "forward") {
		rScat << "forwardInflectionDF<- data.frame(bin,count,variable)\n"; 
		rScat << "forwardStackBarInf <- ggplot(forwardInflectionDF, aes(x = bin, y = count, fill = variable)) + \n"; 
		rScat << "scale_fill_manual(values=c(\"red\", \"blue\", \"chartreuse4\", \"darkgoldenrod1\")) + \n"; 
		rScat << "geom_bar(stat = \"identity\") +\n";
	}
	if (s_direction == "reverse") {
		rScat << "reverseInflectionDF<- data.frame(bin,count,variable)\n"; 
		rScat << "reverseStackBarInf <- ggplot(reverseInflectionDF, aes(x = bin, y = count, fill = variable)) + \n"; 
		rScat << "scale_fill_manual(values=c(\"blue\", \"red\", \"chartreuse4\", \"darkgoldenrod1\")) + \n"; 
		rScat << "geom_bar(stat = \"identity\") +\n";
	}
	rScat << "	ggtitle(\"\") +\n"; 
	rScat << "	labs(x = \"TTI (time)\",size = 10) +\n"; 
	rScat << "	scale_x_continuous(expand = c(0,0), limits = c(0, " << turnoverTimes[i_bamFiles-1] << ")) +\n"; 
	rScat << "	scale_y_continuous(limits = c(0,NA), expand = c(0, 0)) +\n"; 
	rScat << "	theme(\n"; 
	rScat << "		axis.text = element_text(size = 8),\n"; 
	rScat << "		legend.key = element_rect(fill = \"white\"),\n"; 
	rScat << "		legend.background = element_rect(fill = \"white\"),\n"; 
	rScat << "		legend.position=\"none\",\n"; 
	rScat << "		panel.grid.major = element_line(colour = \"white\"),\n"; 
	rScat << "		panel.grid.minor = element_blank(),\n"; 
	rScat << "		panel.background = element_rect(fill = \"white\"),\n"; 
	rScat << "		axis.line.x = element_line(colour = \"black\"),\n"; 
	rScat << "		axis.line.y = element_line(colour = \"black\")\n"; 
	rScat << "	)\n"; 


	// 4 rows, 25 columns, initialized to 0. 
	// rows are: 0) rise/fall 1) hill 2) valley 3) undefined
	int i_hcStackBarArr[4][25]={0};



	for (int i = 1; i < i_peakNumber+1; i++) {	
		if (s_direction == "forward") {
			std::string s = dataArray[i][2*i_bamFiles + 5]; 	 // forward hills coefficient
			if ( (s != "-") &&  (s != "NaN") ) {	 		 // any rise
				double d = stod(dataArray[i][2*i_bamFiles + 5]);
				std::string s_t = dataArray[i][2*i_bamFiles + 3];
			
				if ( (d < 0) && (d >= (d_hcBin*1)) && (s_t == "rise") )	 		   { i_hcStackBarArr[0][0]++; }
				if ( (d < 0) && (d >= (d_hcBin*1)) && (s_t == "hill") )	 		   { i_hcStackBarArr[1][0]++; }
				if ( (d < 0) && (d >= (d_hcBin*1)) && (s_t == "valley") )		  	   { i_hcStackBarArr[2][0]++; }
				if ( (d < 0) && (d >= (d_hcBin*1)) && (s_t == "undefined") )		   { i_hcStackBarArr[3][0]++; }

				if ( (d < (d_hcBin*1)) && (d >= (d_hcBin*2)) && (s_t == "rise") )	  	   { i_hcStackBarArr[0][1]++; }
				if ( (d < (d_hcBin*1)) && (d >= (d_hcBin*2)) && (s_t == "hill") )		   { i_hcStackBarArr[1][1]++; }
				if ( (d < (d_hcBin*1)) && (d >= (d_hcBin*2)) && (s_t == "valley") )	   { i_hcStackBarArr[2][1]++; }
				if ( (d < (d_hcBin*1)) && (d >= (d_hcBin*2)) && (s_t == "undefined") )     { i_hcStackBarArr[3][1]++; }

				if ( (d < (d_hcBin*2)) && (d >= (d_hcBin*3)) && (s_t == "rise") )	  	   { i_hcStackBarArr[0][2]++; }
				if ( (d < (d_hcBin*2)) && (d >= (d_hcBin*3)) && (s_t == "hill") )		   { i_hcStackBarArr[1][2]++; }
				if ( (d < (d_hcBin*2)) && (d >= (d_hcBin*3)) && (s_t == "valley") )	   { i_hcStackBarArr[2][2]++; }
				if ( (d < (d_hcBin*2)) && (d >= (d_hcBin*3)) && (s_t == "undefined") )     { i_hcStackBarArr[3][2]++; }

				if ( (d < (d_hcBin*3)) && (d >= (d_hcBin*4)) && (s_t == "rise") )	   	   { i_hcStackBarArr[0][3]++; }
				if ( (d < (d_hcBin*3)) && (d >= (d_hcBin*4)) && (s_t == "hill") )	  	   { i_hcStackBarArr[1][3]++; }
				if ( (d < (d_hcBin*3)) && (d >= (d_hcBin*4)) && (s_t == "valley") )	   { i_hcStackBarArr[2][3]++; }
				if ( (d < (d_hcBin*3)) && (d >= (d_hcBin*4)) && (s_t == "undefined") )     { i_hcStackBarArr[3][3]++; }

				if ( (d < (d_hcBin*4)) && (d >= (d_hcBin*5)) && (s_t == "rise") )	  	   { i_hcStackBarArr[0][4]++; }
				if ( (d < (d_hcBin*4)) && (d >= (d_hcBin*5)) && (s_t == "hill") )	  	   { i_hcStackBarArr[1][4]++; }
				if ( (d < (d_hcBin*4)) && (d >= (d_hcBin*5)) && (s_t == "valley") )	   { i_hcStackBarArr[2][4]++; }
				if ( (d < (d_hcBin*4)) && (d >= (d_hcBin*5)) && (s_t == "undefined") )     { i_hcStackBarArr[3][4]++; }

				if ( (d < (d_hcBin*5)) && (d >= (d_hcBin*6)) && (s_t == "rise") )	 	   { i_hcStackBarArr[0][5]++; }
				if ( (d < (d_hcBin*5)) && (d >= (d_hcBin*6)) && (s_t == "hill") )	 	   { i_hcStackBarArr[1][5]++; }
				if ( (d < (d_hcBin*5)) && (d >= (d_hcBin*6)) && (s_t == "valley") )	   { i_hcStackBarArr[2][5]++; }
				if ( (d < (d_hcBin*5)) && (d >= (d_hcBin*6)) && (s_t == "undefined") )     { i_hcStackBarArr[3][5]++; }

				if ( (d < (d_hcBin*6)) && (d >= (d_hcBin*7)) && (s_t == "rise") )	   	   { i_hcStackBarArr[0][6]++; }
				if ( (d < (d_hcBin*6)) && (d >= (d_hcBin*7)) && (s_t == "hill") )	 	   { i_hcStackBarArr[1][6]++; }
				if ( (d < (d_hcBin*6)) && (d >= (d_hcBin*7)) && (s_t == "valley") )	   { i_hcStackBarArr[2][6]++; }
				if ( (d < (d_hcBin*6)) && (d >= (d_hcBin*7)) && (s_t == "undefined") )     { i_hcStackBarArr[3][6]++; }

				if ( (d < (d_hcBin*7)) && (d >= (d_hcBin*8)) && (s_t == "rise") )	         { i_hcStackBarArr[0][7]++; }
				if ( (d < (d_hcBin*7)) && (d >= (d_hcBin*8)) && (s_t == "hill") )	         { i_hcStackBarArr[1][7]++; }
				if ( (d < (d_hcBin*7)) && (d >= (d_hcBin*8)) && (s_t == "valley") )	   { i_hcStackBarArr[2][7]++; }
				if ( (d < (d_hcBin*7)) && (d >= (d_hcBin*8)) && (s_t == "undefined") )     { i_hcStackBarArr[3][7]++; }

				if ( (d < (d_hcBin*8)) && (d >= (d_hcBin*9)) && (s_t == "rise") )	         { i_hcStackBarArr[0][8]++; }
				if ( (d < (d_hcBin*8)) && (d >= (d_hcBin*9)) && (s_t == "hill") )	   	   { i_hcStackBarArr[1][8]++; }
				if ( (d < (d_hcBin*8)) && (d >= (d_hcBin*9)) && (s_t == "valley") )	   { i_hcStackBarArr[2][8]++; }
				if ( (d < (d_hcBin*8)) && (d >= (d_hcBin*9)) && (s_t == "undefined") )     { i_hcStackBarArr[3][8]++; }

				if ( (d < (d_hcBin*9)) && (d >= (d_hcBin*10)) && (s_t == "rise") )	   { i_hcStackBarArr[0][9]++; }
				if ( (d < (d_hcBin*9)) && (d >= (d_hcBin*10)) && (s_t == "hill") )	   { i_hcStackBarArr[1][9]++; }
				if ( (d < (d_hcBin*9)) && (d >= (d_hcBin*10)) && (s_t == "valley") )	   { i_hcStackBarArr[2][9]++; }
				if ( (d < (d_hcBin*9)) && (d >= (d_hcBin*10)) && (s_t == "undefined") )    { i_hcStackBarArr[3][9]++; }

				if ( (d < (d_hcBin*10)) && (d >= (d_hcBin*11)) && (s_t == "rise") )	   { i_hcStackBarArr[0][10]++; }
				if ( (d < (d_hcBin*10)) && (d >= (d_hcBin*11)) && (s_t == "hill") )	   { i_hcStackBarArr[1][10]++; }
				if ( (d < (d_hcBin*10)) && (d >= (d_hcBin*11)) && (s_t == "valley") )	   { i_hcStackBarArr[2][10]++; }
				if ( (d < (d_hcBin*10)) && (d >= (d_hcBin*11)) && (s_t == "undefined") )   { i_hcStackBarArr[3][10]++; }

				if ( (d < (d_hcBin*11)) && (d >= (d_hcBin*12)) && (s_t == "rise") )	   { i_hcStackBarArr[0][11]++; }
				if ( (d < (d_hcBin*11)) && (d >= (d_hcBin*12)) && (s_t == "hill") )	   { i_hcStackBarArr[1][11]++; }
				if ( (d < (d_hcBin*11)) && (d >= (d_hcBin*12)) && (s_t == "valley") )	   { i_hcStackBarArr[2][11]++; }
				if ( (d < (d_hcBin*11)) && (d >= (d_hcBin*12)) && (s_t == "undefined") )   { i_hcStackBarArr[3][11]++; }

				if ( (d < (d_hcBin*12)) && (d >= (d_hcBin*13)) && (s_t == "rise") )	   { i_hcStackBarArr[0][12]++; }
				if ( (d < (d_hcBin*12)) && (d >= (d_hcBin*13)) && (s_t == "hill") )	   { i_hcStackBarArr[1][12]++; }
				if ( (d < (d_hcBin*12)) && (d >= (d_hcBin*13)) && (s_t == "valley") )	   { i_hcStackBarArr[2][12]++; }
				if ( (d < (d_hcBin*12)) && (d >= (d_hcBin*13)) && (s_t == "undefined") )   { i_hcStackBarArr[3][12]++; }

				if ( (d < (d_hcBin*13)) && (d >= (d_hcBin*14)) && (s_t == "rise") )	   { i_hcStackBarArr[0][13]++; }
				if ( (d < (d_hcBin*13)) && (d >= (d_hcBin*14)) && (s_t == "hill") )	   { i_hcStackBarArr[1][13]++; }
				if ( (d < (d_hcBin*13)) && (d >= (d_hcBin*14)) && (s_t == "valley") )	   { i_hcStackBarArr[2][13]++; }
				if ( (d < (d_hcBin*13)) && (d >= (d_hcBin*14)) && (s_t == "undefined") )   { i_hcStackBarArr[3][13]++; }

				if ( (d < (d_hcBin*14)) && (d >= (d_hcBin*15)) && (s_t == "rise") )	   { i_hcStackBarArr[0][14]++; }
				if ( (d < (d_hcBin*14)) && (d >= (d_hcBin*15)) && (s_t == "hill") )	   { i_hcStackBarArr[1][14]++; }
				if ( (d < (d_hcBin*14)) && (d >= (d_hcBin*15)) && (s_t == "valley") )	   { i_hcStackBarArr[2][14]++; }
				if ( (d < (d_hcBin*14)) && (d >= (d_hcBin*15)) && (s_t == "undefined") )   { i_hcStackBarArr[3][14]++; }

				if ( (d < (d_hcBin*15)) && (d >= (d_hcBin*16)) && (s_t == "rise") )	   { i_hcStackBarArr[0][15]++; }
				if ( (d < (d_hcBin*15)) && (d >= (d_hcBin*16)) && (s_t == "hill") )	   { i_hcStackBarArr[1][15]++; }
				if ( (d < (d_hcBin*15)) && (d >= (d_hcBin*16)) && (s_t == "valley") )	   { i_hcStackBarArr[2][15]++; }
				if ( (d < (d_hcBin*15)) && (d >= (d_hcBin*16)) && (s_t == "undefined") )   { i_hcStackBarArr[3][15]++; }

				if ( (d < (d_hcBin*16)) && (d >= (d_hcBin*17)) && (s_t == "rise") )	   { i_hcStackBarArr[0][16]++; }
				if ( (d < (d_hcBin*16)) && (d >= (d_hcBin*17)) && (s_t == "hill") )	   { i_hcStackBarArr[1][16]++; }
				if ( (d < (d_hcBin*16)) && (d >= (d_hcBin*17)) && (s_t == "valley") )	   { i_hcStackBarArr[2][16]++; }
				if ( (d < (d_hcBin*16)) && (d >= (d_hcBin*17)) && (s_t == "undefined") )   { i_hcStackBarArr[3][16]++; }

				if ( (d < (d_hcBin*17)) && (d >= (d_hcBin*18)) && (s_t == "rise") )	   { i_hcStackBarArr[0][17]++; }
				if ( (d < (d_hcBin*17)) && (d >= (d_hcBin*18)) && (s_t == "hill") )	   { i_hcStackBarArr[1][17]++; }
				if ( (d < (d_hcBin*17)) && (d >= (d_hcBin*18)) && (s_t == "valley") )	   { i_hcStackBarArr[2][17]++; }
				if ( (d < (d_hcBin*17)) && (d >= (d_hcBin*18)) && (s_t == "undefined") )   { i_hcStackBarArr[3][17]++; }

				if ( (d < (d_hcBin*18)) && (d >= (d_hcBin*19)) && (s_t == "rise") )	   { i_hcStackBarArr[0][18]++; }
				if ( (d < (d_hcBin*18)) && (d >= (d_hcBin*19)) && (s_t == "hill") )	   { i_hcStackBarArr[1][18]++; }
				if ( (d < (d_hcBin*18)) && (d >= (d_hcBin*19)) && (s_t == "valley") )	   { i_hcStackBarArr[2][18]++; }
				if ( (d < (d_hcBin*18)) && (d >= (d_hcBin*19)) && (s_t == "undefined") )   { i_hcStackBarArr[3][18]++; }

				if ( (d < (d_hcBin*19)) && (d >= (d_hcBin*20)) && (s_t == "rise") )	   { i_hcStackBarArr[0][19]++; }
				if ( (d < (d_hcBin*19)) && (d >= (d_hcBin*20)) && (s_t == "hill") )	   { i_hcStackBarArr[1][19]++; }
				if ( (d < (d_hcBin*19)) && (d >= (d_hcBin*20)) && (s_t == "valley") )	   { i_hcStackBarArr[2][19]++; }
				if ( (d < (d_hcBin*19)) && (d >= (d_hcBin*20)) && (s_t == "undefined") )   { i_hcStackBarArr[3][19]++; }

				if ( (d < (d_hcBin*20)) && (d >= (d_hcBin*21)) && (s_t == "rise") )	   { i_hcStackBarArr[0][20]++; }
				if ( (d < (d_hcBin*20)) && (d >= (d_hcBin*21)) && (s_t == "hill") )	   { i_hcStackBarArr[1][20]++; }
				if ( (d < (d_hcBin*20)) && (d >= (d_hcBin*21)) && (s_t == "valley") )	   { i_hcStackBarArr[2][20]++; }
				if ( (d < (d_hcBin*20)) && (d >= (d_hcBin*21)) && (s_t == "undefined") )   { i_hcStackBarArr[3][20]++; }

				if ( (d < (d_hcBin*21)) && (d >= (d_hcBin*22)) && (s_t == "rise") )	   { i_hcStackBarArr[0][21]++; }
				if ( (d < (d_hcBin*21)) && (d >= (d_hcBin*22)) && (s_t == "hill") )	   { i_hcStackBarArr[1][21]++; }
				if ( (d < (d_hcBin*21)) && (d >= (d_hcBin*22)) && (s_t == "valley") )	   { i_hcStackBarArr[2][21]++; }
				if ( (d < (d_hcBin*21)) && (d >= (d_hcBin*22)) && (s_t == "undefined") )   { i_hcStackBarArr[3][21]++; }

				if ( (d < (d_hcBin*22)) && (d >= (d_hcBin*23)) && (s_t == "rise") )	   { i_hcStackBarArr[0][22]++; }
				if ( (d < (d_hcBin*22)) && (d >= (d_hcBin*23)) && (s_t == "hill") )	   { i_hcStackBarArr[1][22]++; }
				if ( (d < (d_hcBin*22)) && (d >= (d_hcBin*23)) && (s_t == "valley") )	   { i_hcStackBarArr[2][22]++; }
				if ( (d < (d_hcBin*22)) && (d >= (d_hcBin*23)) && (s_t == "undefined") )   { i_hcStackBarArr[3][22]++; }

				if ( (d < (d_hcBin*23)) && (d >= (d_hcBin*24)) && (s_t == "rise") )	   { i_hcStackBarArr[0][23]++; }
				if ( (d < (d_hcBin*23)) && (d >= (d_hcBin*24)) && (s_t == "hill") )	   { i_hcStackBarArr[1][23]++; }
				if ( (d < (d_hcBin*23)) && (d >= (d_hcBin*24)) && (s_t == "valley") )	   { i_hcStackBarArr[2][23]++; }
				if ( (d < (d_hcBin*23)) && (d >= (d_hcBin*24)) && (s_t == "undefined") )   { i_hcStackBarArr[3][23]++; }

				if ( (d < (d_hcBin*24)) && (d >= (d_hcBin*25)) && (s_t == "rise") )	   { i_hcStackBarArr[0][24]++; }
				if ( (d < (d_hcBin*24)) && (d >= (d_hcBin*25)) && (s_t == "hill") )	   { i_hcStackBarArr[1][24]++; }
				if ( (d < (d_hcBin*24)) && (d >= (d_hcBin*25)) && (s_t == "valley") )	   { i_hcStackBarArr[2][24]++; }
				if ( (d < (d_hcBin*24)) && (d >= (d_hcBin*25)) && (s_t == "undefined") )   { i_hcStackBarArr[3][24]++; }
			}
		}
		if (s_direction == "reverse") {
			std::string s = dataArray[i][2*i_bamFiles + 9]; 	 // reverse hills coefficient
			if ( (s != "-") &&  (s != "NaN") ) {	 		 // any fall
				double d = stod(dataArray[i][2*i_bamFiles + 9]); 
				std::string s_t = dataArray[i][2*i_bamFiles + 3];
			
				if ( (d > 0) && (d <= (d_hcBin*1)) && (s_t == "fall") )	 		   { i_hcStackBarArr[0][0]++; }
				if ( (d > 0) && (d <= (d_hcBin*1)) && (s_t == "hill") )	 		   { i_hcStackBarArr[1][0]++; }
				if ( (d > 0) && (d <= (d_hcBin*1)) && (s_t == "valley") )		 	   { i_hcStackBarArr[2][0]++; }
				if ( (d > 0) && (d <= (d_hcBin*1)) && (s_t == "undefined") )		   { i_hcStackBarArr[3][0]++; }

				if ( (d > (d_hcBin*1)) && (d <= (d_hcBin*2)) && (s_t == "fall") )	  	   { i_hcStackBarArr[0][1]++; }
				if ( (d > (d_hcBin*1)) && (d <= (d_hcBin*2)) && (s_t == "hill") )	  	   { i_hcStackBarArr[1][1]++; }
				if ( (d > (d_hcBin*1)) && (d <= (d_hcBin*2)) && (s_t == "valley") )	   { i_hcStackBarArr[2][1]++; }
				if ( (d > (d_hcBin*1)) && (d <= (d_hcBin*2)) && (s_t == "undefined") )     { i_hcStackBarArr[3][1]++; }

				if ( (d > (d_hcBin*2)) && (d <= (d_hcBin*3)) && (s_t == "fall") )	  	   { i_hcStackBarArr[0][2]++; }
				if ( (d > (d_hcBin*2)) && (d <= (d_hcBin*3)) && (s_t == "hill") )	   	   { i_hcStackBarArr[1][2]++; }
				if ( (d > (d_hcBin*2)) && (d <= (d_hcBin*3)) && (s_t == "valley") )	   { i_hcStackBarArr[2][2]++; }
				if ( (d > (d_hcBin*2)) && (d <= (d_hcBin*3)) && (s_t == "undefined") )     { i_hcStackBarArr[3][2]++; }

				if ( (d > (d_hcBin*3)) && (d <= (d_hcBin*4)) && (s_t == "fall") )	 	   { i_hcStackBarArr[0][3]++; }
				if ( (d > (d_hcBin*3)) && (d <= (d_hcBin*4)) && (s_t == "hill") )	  	   { i_hcStackBarArr[1][3]++; }
				if ( (d > (d_hcBin*3)) && (d <= (d_hcBin*4)) && (s_t == "valley") )	   { i_hcStackBarArr[2][3]++; }
				if ( (d > (d_hcBin*3)) && (d <= (d_hcBin*4)) && (s_t == "undefined") )     { i_hcStackBarArr[3][3]++; }

				if ( (d > (d_hcBin*4)) && (d <= (d_hcBin*5)) && (s_t == "fall") )	 	   { i_hcStackBarArr[0][4]++; }
				if ( (d > (d_hcBin*4)) && (d <= (d_hcBin*5)) && (s_t == "hill") )	 	   { i_hcStackBarArr[1][4]++; }
				if ( (d > (d_hcBin*4)) && (d <= (d_hcBin*5)) && (s_t == "valley") )	   { i_hcStackBarArr[2][4]++; }
				if ( (d > (d_hcBin*4)) && (d <= (d_hcBin*5)) && (s_t == "undefined") )     { i_hcStackBarArr[3][4]++; }

				if ( (d > (d_hcBin*5)) && (d <= (d_hcBin*6)) && (s_t == "fall") )	  	   { i_hcStackBarArr[0][5]++; }
				if ( (d > (d_hcBin*5)) && (d <= (d_hcBin*6)) && (s_t == "hill") )	  	   { i_hcStackBarArr[1][5]++; }
				if ( (d > (d_hcBin*5)) && (d <= (d_hcBin*6)) && (s_t == "valley") )	   { i_hcStackBarArr[2][5]++; }
				if ( (d > (d_hcBin*5)) && (d <= (d_hcBin*6)) && (s_t == "undefined") )     { i_hcStackBarArr[3][5]++; }

				if ( (d > (d_hcBin*6)) && (d <= (d_hcBin*7)) && (s_t == "fall") )	 	   { i_hcStackBarArr[0][6]++; }
				if ( (d > (d_hcBin*6)) && (d <= (d_hcBin*7)) && (s_t == "hill") )	   	   { i_hcStackBarArr[1][6]++; }
				if ( (d > (d_hcBin*6)) && (d <= (d_hcBin*7)) && (s_t == "valley") )	   { i_hcStackBarArr[2][6]++; }
				if ( (d > (d_hcBin*6)) && (d <= (d_hcBin*7)) && (s_t == "undefined") )     { i_hcStackBarArr[3][6]++; }

				if ( (d > (d_hcBin*7)) && (d <= (d_hcBin*8)) && (s_t == "fall") )	   	   { i_hcStackBarArr[0][7]++; }
				if ( (d > (d_hcBin*7)) && (d <= (d_hcBin*8)) && (s_t == "hill") )	   	   { i_hcStackBarArr[1][7]++; }
				if ( (d > (d_hcBin*7)) && (d <= (d_hcBin*8)) && (s_t == "valley") )	   { i_hcStackBarArr[2][7]++; }
				if ( (d > (d_hcBin*7)) && (d <= (d_hcBin*8)) && (s_t == "undefined") )     { i_hcStackBarArr[3][7]++; }

				if ( (d > (d_hcBin*8)) && (d <= (d_hcBin*9)) && (s_t == "fall") )	   	   { i_hcStackBarArr[0][8]++; }
				if ( (d > (d_hcBin*8)) && (d <= (d_hcBin*9)) && (s_t == "hill") )	  	   { i_hcStackBarArr[1][8]++; }
				if ( (d > (d_hcBin*8)) && (d <= (d_hcBin*9)) && (s_t == "valley") )	   { i_hcStackBarArr[2][8]++; }
				if ( (d > (d_hcBin*8)) && (d <= (d_hcBin*9)) && (s_t == "undefined") )     { i_hcStackBarArr[3][8]++; }

				if ( (d > (d_hcBin*9)) && (d <= (d_hcBin*10)) && (s_t == "fall") )	   { i_hcStackBarArr[0][9]++; }
				if ( (d > (d_hcBin*9)) && (d <= (d_hcBin*10)) && (s_t == "hill") )	   { i_hcStackBarArr[1][9]++; }
				if ( (d > (d_hcBin*9)) && (d <= (d_hcBin*10)) && (s_t == "valley") )	   { i_hcStackBarArr[2][9]++; }
				if ( (d > (d_hcBin*9)) && (d <= (d_hcBin*10)) && (s_t == "undefined") )    { i_hcStackBarArr[3][9]++; }

				if ( (d > (d_hcBin*10)) && (d <= (d_hcBin*11)) && (s_t == "fall") )	   { i_hcStackBarArr[0][10]++; }
				if ( (d > (d_hcBin*10)) && (d <= (d_hcBin*11)) && (s_t == "hill") )	   { i_hcStackBarArr[1][10]++; }
				if ( (d > (d_hcBin*10)) && (d <= (d_hcBin*11)) && (s_t == "valley") )	   { i_hcStackBarArr[2][10]++; }
				if ( (d > (d_hcBin*10)) && (d <= (d_hcBin*11)) && (s_t == "undefined") )   { i_hcStackBarArr[3][10]++; }

				if ( (d > (d_hcBin*11)) && (d <= (d_hcBin*12)) && (s_t == "fall") )	   { i_hcStackBarArr[0][11]++; }
				if ( (d > (d_hcBin*11)) && (d <= (d_hcBin*12)) && (s_t == "hill") )	   { i_hcStackBarArr[1][11]++; }
				if ( (d > (d_hcBin*11)) && (d <= (d_hcBin*12)) && (s_t == "valley") )	   { i_hcStackBarArr[2][11]++; }
				if ( (d > (d_hcBin*11)) && (d <= (d_hcBin*12)) && (s_t == "undefined") )   { i_hcStackBarArr[3][11]++; }

				if ( (d > (d_hcBin*12)) && (d <= (d_hcBin*13)) && (s_t == "fall") )	   { i_hcStackBarArr[0][12]++; }
				if ( (d > (d_hcBin*12)) && (d <= (d_hcBin*13)) && (s_t == "hill") )	   { i_hcStackBarArr[1][12]++; }
				if ( (d > (d_hcBin*12)) && (d <= (d_hcBin*13)) && (s_t == "valley") )	   { i_hcStackBarArr[2][12]++; }
				if ( (d > (d_hcBin*12)) && (d <= (d_hcBin*13)) && (s_t == "undefined") )   { i_hcStackBarArr[3][12]++; }

				if ( (d > (d_hcBin*13)) && (d <= (d_hcBin*14)) && (s_t == "fall") )	   { i_hcStackBarArr[0][13]++; }
				if ( (d > (d_hcBin*13)) && (d <= (d_hcBin*14)) && (s_t == "hill") )	   { i_hcStackBarArr[1][13]++; }
				if ( (d > (d_hcBin*13)) && (d <= (d_hcBin*14)) && (s_t == "valley") )	   { i_hcStackBarArr[2][13]++; }
				if ( (d > (d_hcBin*13)) && (d <= (d_hcBin*14)) && (s_t == "undefined") )   { i_hcStackBarArr[3][13]++; }

				if ( (d > (d_hcBin*14)) && (d <= (d_hcBin*15)) && (s_t == "fall") )	   { i_hcStackBarArr[0][14]++; }
				if ( (d > (d_hcBin*14)) && (d <= (d_hcBin*15)) && (s_t == "hill") )	   { i_hcStackBarArr[1][14]++; }
				if ( (d > (d_hcBin*14)) && (d <= (d_hcBin*15)) && (s_t == "valley") )	   { i_hcStackBarArr[2][14]++; }
				if ( (d > (d_hcBin*14)) && (d <= (d_hcBin*15)) && (s_t == "undefined") )   { i_hcStackBarArr[3][14]++; }

				if ( (d > (d_hcBin*15)) && (d <= (d_hcBin*16)) && (s_t == "fall") )	   { i_hcStackBarArr[0][15]++; }
				if ( (d > (d_hcBin*15)) && (d <= (d_hcBin*16)) && (s_t == "hill") )	   { i_hcStackBarArr[1][15]++; }
				if ( (d > (d_hcBin*15)) && (d <= (d_hcBin*16)) && (s_t == "valley") )	   { i_hcStackBarArr[2][15]++; }
				if ( (d > (d_hcBin*15)) && (d <= (d_hcBin*16)) && (s_t == "undefined") )   { i_hcStackBarArr[3][15]++; }

				if ( (d > (d_hcBin*16)) && (d <= (d_hcBin*17)) && (s_t == "fall") )	   { i_hcStackBarArr[0][16]++; }
				if ( (d > (d_hcBin*16)) && (d <= (d_hcBin*17)) && (s_t == "hill") )	   { i_hcStackBarArr[1][16]++; }
				if ( (d > (d_hcBin*16)) && (d <= (d_hcBin*17)) && (s_t == "valley") )	   { i_hcStackBarArr[2][16]++; }
				if ( (d > (d_hcBin*16)) && (d <= (d_hcBin*17)) && (s_t == "undefined") )   { i_hcStackBarArr[3][16]++; }

				if ( (d > (d_hcBin*17)) && (d <= (d_hcBin*18)) && (s_t == "fall") )	   { i_hcStackBarArr[0][17]++; }
				if ( (d > (d_hcBin*17)) && (d <= (d_hcBin*18)) && (s_t == "hill") )	   { i_hcStackBarArr[1][17]++; }
				if ( (d > (d_hcBin*17)) && (d <= (d_hcBin*18)) && (s_t == "valley") )	   { i_hcStackBarArr[2][17]++; }
				if ( (d > (d_hcBin*17)) && (d <= (d_hcBin*18)) && (s_t == "undefined") )   { i_hcStackBarArr[3][17]++; }

				if ( (d > (d_hcBin*18)) && (d <= (d_hcBin*19)) && (s_t == "fall") )	   { i_hcStackBarArr[0][18]++; }
				if ( (d > (d_hcBin*18)) && (d <= (d_hcBin*19)) && (s_t == "hill") )	   { i_hcStackBarArr[1][18]++; }
				if ( (d > (d_hcBin*18)) && (d <= (d_hcBin*19)) && (s_t == "valley") )	   { i_hcStackBarArr[2][18]++; }
				if ( (d > (d_hcBin*18)) && (d <= (d_hcBin*19)) && (s_t == "undefined") )   { i_hcStackBarArr[3][18]++; }

				if ( (d > (d_hcBin*19)) && (d <= (d_hcBin*20)) && (s_t == "fall") )	   { i_hcStackBarArr[0][19]++; }
				if ( (d > (d_hcBin*19)) && (d <= (d_hcBin*20)) && (s_t == "hill") )	   { i_hcStackBarArr[1][19]++; }
				if ( (d > (d_hcBin*19)) && (d <= (d_hcBin*20)) && (s_t == "valley") )	   { i_hcStackBarArr[2][19]++; }
				if ( (d > (d_hcBin*19)) && (d <= (d_hcBin*20)) && (s_t == "undefined") )   { i_hcStackBarArr[3][19]++; }

				if ( (d > (d_hcBin*20)) && (d <= (d_hcBin*21)) && (s_t == "fall") )	   { i_hcStackBarArr[0][20]++; }
				if ( (d > (d_hcBin*20)) && (d <= (d_hcBin*21)) && (s_t == "hill") )	   { i_hcStackBarArr[1][20]++; }
				if ( (d > (d_hcBin*20)) && (d <= (d_hcBin*21)) && (s_t == "valley") )	   { i_hcStackBarArr[2][20]++; }
				if ( (d > (d_hcBin*20)) && (d <= (d_hcBin*21)) && (s_t == "undefined") )   { i_hcStackBarArr[3][20]++; }

				if ( (d > (d_hcBin*21)) && (d <= (d_hcBin*22)) && (s_t == "fall") )	   { i_hcStackBarArr[0][21]++; }
				if ( (d > (d_hcBin*21)) && (d <= (d_hcBin*22)) && (s_t == "hill") )	   { i_hcStackBarArr[1][21]++; }
				if ( (d > (d_hcBin*21)) && (d <= (d_hcBin*22)) && (s_t == "valley") )	   { i_hcStackBarArr[2][21]++; }
				if ( (d > (d_hcBin*21)) && (d <= (d_hcBin*22)) && (s_t == "undefined") )   { i_hcStackBarArr[3][21]++; }

				if ( (d > (d_hcBin*22)) && (d <= (d_hcBin*23)) && (s_t == "fall") )	   { i_hcStackBarArr[0][22]++; }
				if ( (d > (d_hcBin*22)) && (d <= (d_hcBin*23)) && (s_t == "hill") )	   { i_hcStackBarArr[1][22]++; }
				if ( (d > (d_hcBin*22)) && (d <= (d_hcBin*23)) && (s_t == "valley") )	   { i_hcStackBarArr[2][22]++; }
				if ( (d > (d_hcBin*22)) && (d <= (d_hcBin*23)) && (s_t == "undefined") )   { i_hcStackBarArr[3][22]++; }

				if ( (d > (d_hcBin*23)) && (d <= (d_hcBin*24)) && (s_t == "fall") )	   { i_hcStackBarArr[0][23]++; }
				if ( (d > (d_hcBin*23)) && (d <= (d_hcBin*24)) && (s_t == "hill") )	   { i_hcStackBarArr[1][23]++; }
				if ( (d > (d_hcBin*23)) && (d <= (d_hcBin*24)) && (s_t == "valley") )	   { i_hcStackBarArr[2][23]++; }
				if ( (d > (d_hcBin*23)) && (d <= (d_hcBin*24)) && (s_t == "undefined") )   { i_hcStackBarArr[3][23]++; }

				if ( (d > (d_hcBin*24)) && (d <= (d_hcBin*25)) && (s_t == "fall") )	   { i_hcStackBarArr[0][24]++; }
				if ( (d > (d_hcBin*24)) && (d <= (d_hcBin*25)) && (s_t == "hill") )	   { i_hcStackBarArr[1][24]++; }
				if ( (d > (d_hcBin*24)) && (d <= (d_hcBin*25)) && (s_t == "valley") )	   { i_hcStackBarArr[2][24]++; }
				if ( (d > (d_hcBin*24)) && (d <= (d_hcBin*25)) && (s_t == "undefined") )   { i_hcStackBarArr[3][24]++; }
			}
		}
	}









	// rows are: 0) rise/fall 1) hill 2) valley 3) undefined
	rScat << "hcBin <- c(";
	for (int j = 0; j < 3; j++) {	
		for (int i = 0; i < 25; i++) {	
			rScat << (d_hcBin*i) + d_hcBin/2 << ", ";
		}
	}
	for (int i = 0; i < 24; i++) {	
		rScat << (d_hcBin*i) + d_hcBin/2 << ", ";
	}
	rScat << (d_hcBin*24) + d_hcBin/2 << ")\n";

	// print count
	rScat << "count <- c(";
	for (int i = 0; i < 3; i++) {	
		for (int j = 0; j < 25; j++) {	
			rScat << i_hcStackBarArr[i][j] << ", ";
		}
	}
	for (int i = 0; i < 24; i++) {	
		rScat << i_hcStackBarArr[3][i] << ", ";
	}
	rScat << i_hcStackBarArr[3][24] << ")\n";
	
	// print variable 
	rScat << "variable <- c(\n";
	if (s_direction == "forward") {
		for (int i = 0; i < 25; i++) {	
			rScat << "\"rise\", ";
		}
		for (int i = 0; i < 25; i++) {	
			rScat << "\"hill\", ";
		}
		for (int i = 0; i < 25; i++) {	
			rScat << "\"valley\", ";
		}
		for (int i = 0; i < 24; i++) {	
			rScat << "\"undefined rise\", ";
		}
		rScat << "\"undefined rise\")\n";
	}
	if (s_direction == "reverse") {
		for (int i = 0; i < 25; i++) {	
			rScat << "\"fall\", ";
		}
		for (int i = 0; i < 25; i++) {	
			rScat << "\"hill\", ";
		}
		for (int i = 0; i < 25; i++) {	
			rScat << "\"valley\", ";
		}
		for (int i = 0; i < 24; i++) {	
			rScat << "\"undefined fall\", ";
		}
		rScat << "\"undefined fall\")\n";
	}

	if (s_direction == "forward") {
		rScat << "forwardHcDF<- data.frame(hcBin,count,variable)\n"; 
		rScat << "forwardStackBarHC <- ggplot(forwardHcDF, aes(x = hcBin, y = count, fill = variable)) + \n"; 
		rScat << "scale_fill_manual(values=c(\"red\", \"blue\", \"chartreuse4\", \"darkgoldenrod1\")) + \n"; 
		rScat << "coord_flip() + \n"; 
		rScat << "geom_bar(stat = \"identity\") +\n";
	}
	if (s_direction == "reverse") {
		rScat << "reverseHcDF<- data.frame(hcBin,count,variable)\n"; 
		rScat << "reverseStackBarHC <- ggplot(reverseHcDF, aes(x = hcBin, y = count, fill = variable)) + \n"; 
		rScat << "scale_fill_manual(values=c(\"blue\", \"red\", \"chartreuse4\", \"darkgoldenrod1\")) + \n"; 
		rScat << "coord_flip() + \n"; 
		rScat << "geom_bar(stat = \"identity\") +\n";
	}
	rScat << "	labs(x = \"Hills Coefficient\",size = 10) +\n"; 
	if (s_direction == "forward") {
		if (d_hcVec.size() < 10) { rScat << "	scale_x_continuous(expand = c(0,0), limits = c(-1, 0)) +\n"; }
		else { rScat << "	scale_x_continuous(expand = c(0,0), limits = c(" << d_hcVec[i_index10per] << ", 0)) +\n"; }

	}
	if (s_direction == "reverse") {
		if (d_hcVec.size() < 10) { rScat << "	scale_x_continuous(expand = c(0,0), limits = c(0, 1)) +\n"; } 
		else { rScat << "	scale_x_continuous(expand = c(0,0), limits = c(0, " << d_hcVec[i_index90per] << ")) +\n"; }  
	}
	rScat << "	ggtitle(\"\") +\n"; 
	rScat << "	scale_y_continuous(limits = c(0,NA), expand = c(0, 0)) +\n"; 
	rScat << "	theme(\n"; 
	rScat << "		axis.text = element_text(size = 8),\n"; 
	rScat << "		legend.key = element_rect(fill = \"white\"),\n"; 
	rScat << "		legend.background = element_rect(fill = \"white\"),\n"; 
	rScat << "		legend.position=\"none\",\n"; 
	rScat << "		panel.grid.major = element_line(colour = \"white\"),\n"; 
	rScat << "		panel.grid.minor = element_blank(),\n"; 
	rScat << "		panel.background = element_rect(fill = \"white\"),\n"; 
	rScat << "		axis.line.x = element_line(colour = \"black\"),\n"; 
	rScat << "		axis.line.y = element_line(colour = \"black\")\n"; 
	rScat << "	)\n"; 



	if (b_model) {
		if (s_direction == "reverse") {
			rScat << "d = data.frame(x1=c(1,1,1,1), x2=c(2,2,2,2), y1=c(1,2.5,4,5.5), y2=c(2,3.5,5,6.5), t=c('undefined fall','valley','hill','fall'))\n"; 
			rScat << "legRev <- ggplot() + \n"; 
			rScat << "scale_y_continuous(expand = c(0,0), limits = c(-1, 7.5)) + \n"; 
			rScat << "scale_x_continuous(expand = c(0,0), limits = c(-1, 7)) + \n"; 
			rScat << "geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color=\"black\") +\n"; 
			rScat << "scale_fill_manual(values=c(\"blue\", \"red\", \"chartreuse4\", \"darkgoldenrod1\")) +\n"; 
			rScat << "geom_text(data=d, aes(x=2.5, y=y1+(y2-y1)/2, label=t), size=3, hjust=0) +\n"; 
			rScat << "	theme( \n"; 
			rScat << "		axis.text = element_text(size = 8, colour = \"white\"),\n"; 
			rScat << "		legend.key = element_rect(fill = \"white\"), \n"; 
			rScat << "		legend.background = element_rect(fill = \"white\"),\n"; 
			rScat << "		panel.grid.major = element_line(colour = \"white\"), \n"; 
			rScat << "		legend.position = \"\", \n"; 
			rScat << "		panel.grid.minor = element_blank(),\n"; 
			rScat << "		panel.background = element_rect(fill = \"white\"),\n"; 
			rScat << "		axis.line.x = element_line(colour = \"white\"),\n"; 
			rScat << "		axis.line.y = element_line(colour = \"white\"),\n"; 
			rScat << "		axis.title.x=element_blank(),\n"; 
			rScat << "		axis.text.x=element_blank(),\n"; 
			rScat << "		axis.ticks.x=element_blank(),\n"; 
			rScat << "		axis.title.y=element_blank(),\n"; 
			rScat << "		axis.text.y=element_blank(),\n"; 
			rScat << "		axis.ticks.y=element_blank()\n"; 
			rScat << "		) \n"; 
		}
		if (s_direction == "forward") {
			rScat << "d = data.frame(x1=c(1,1,1,1), x2=c(2,2,2,2), y1=c(1,2.5,4,5.5), y2=c(2,3.5,5,6.5), t=c('undefined rise','valley','hill','rise'))\n"; 
			rScat << "legFor <- ggplot() + \n"; 
			rScat << "scale_y_continuous(expand = c(0,0), limits = c(-1, 7.5)) + \n"; 
			rScat << "scale_x_continuous(expand = c(0,0), limits = c(-1, 7)) + \n"; 
			rScat << "geom_rect(data=d, mapping=aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=t), color=\"black\") +\n"; 
			rScat << "scale_fill_manual(values=c(\"red\", \"blue\", \"chartreuse4\", \"darkgoldenrod1\")) +\n"; 
			rScat << "geom_text(data=d, aes(x=2.5, y=y1+(y2-y1)/2, label=t), size=3, hjust=0) +\n"; 
			rScat << "	theme( \n"; 
			rScat << "		axis.text = element_text(size = 8, colour = \"white\"),\n"; 
			rScat << "		legend.key = element_rect(fill = \"white\"), \n"; 
			rScat << "		legend.background = element_rect(fill = \"white\"),\n"; 
			rScat << "		panel.grid.major = element_line(colour = \"white\"), \n"; 
			rScat << "		legend.position = \"\", \n"; 
			rScat << "		panel.grid.minor = element_blank(),\n"; 
			rScat << "		panel.background = element_rect(fill = \"white\"),\n"; 
			rScat << "		axis.line.x = element_line(colour = \"white\"),\n"; 
			rScat << "		axis.line.y = element_line(colour = \"white\"),\n"; 
			rScat << "		axis.title.x=element_blank(),\n"; 
			rScat << "		axis.text.x=element_blank(),\n"; 
			rScat << "		axis.ticks.x=element_blank(),\n"; 
			rScat << "		axis.title.y=element_blank(),\n"; 
			rScat << "		axis.text.y=element_blank(),\n"; 
			rScat << "		axis.ticks.y=element_blank()\n"; 
			rScat << "		) \n"; 
		}
	}

	d_hcVec.clear();
	rScat.close();
	csvScatInf.close();
	csvScatHC.close();
	csvScatID.close();

}
