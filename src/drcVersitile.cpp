#include <fstream>
#include <iostream>
#include <string>
#include <regex>
#include <vector>
using namespace std;

std::string drcR(int i_loopNumber, std::string s_peakExtention, bool b_L4flag);

void drcVersitile(int i_drcPar, int i_maxThreads, int i_loopNumber, int i_bamFiles, int turnoverTimes[], std::vector<double> &turnoverArr, std::vector<bool> &b_typeArr, int i_minIndex, int i_maxIndex, double d_satthresh, std::vector<bool> &b_forwardArr, std::vector<bool> &b_reverseArr, std::vector<bool> &b_peak1Arr, std::vector<bool> &b_peak2Arr, std::vector<bool> &b_dip1Arr, std::vector<bool> &b_dip2Arr, bool b_model, bool b_L4flag, std::string s_name) 
{

	// write R script
	ofstream drcInflection;
	//drcInflection.open ("drcVersitile." + s_name + "R", std::ios::app);

	if (i_maxThreads == 1) { drcInflection.open ("drcVersitile." + s_name + ".R", std::ios::app); }		//XXX NEW
	else {																//XXX NEW
		if (i_loopNumber < i_drcPar*(i_maxThreads-1)) {									//XXX NEW
			for (int i = 1; i < i_maxThreads; i++) {									//XXX NEW
				if (i_loopNumber < i_drcPar*i) {									//XXX NEW
					string s_drcRfileName = "drcVersitile." + s_name + to_string(i) + ".R";		//XXX NEW
					drcInflection.open (s_drcRfileName, std::ios::app);					//XXX NEW
					break;												//XXX NEW
				}														//XXX NEW
			}															//XXX NEW
		} else { 															//XXX NEW
			string s_drcRfileName = "drcVersitile." + s_name + to_string(i_maxThreads) + ".R";		//XXX NEW
			drcInflection.open (s_drcRfileName, std::ios::app);							//XXX NEW
		}																//XXX NEW
																		//XXX NEW
	}																	//XXX NEW


	int i_start1 = 0;
	int i_end1 = 0;
	int i_start2 = 0;
	int i_end2 = 0;
	
	bool b_allSame = true;
	double d_samCheck = turnoverArr[0];
	//check if all depth values are zero (R can't handle this
	for (int i = 1; i < i_bamFiles; i++) {
		if (d_samCheck != turnoverArr[i]) { b_allSame = false; break;} 	
	}

	if (b_allSame) {
		drcInflection << "print(paste(\"peak" << i_loopNumber << "\"));\n";
		drcInflection << "print(paste(\"e:(Intercept) ERROR\")); print(paste(\"d:(Intercept) ERROR\")); print(paste(\"c:(Intercept) ERROR\"));" <<
			" print(paste(\"b:(Intercept) ERROR\")); print(paste(\"f:(Intercept) ERROR\")); print(paste(\" NaN (-1 degrees of freedom)\"))\n";
	} else {
		if ( ((b_typeArr == b_forwardArr) || (b_typeArr == b_reverseArr)) && (b_model == true) ) {	//forward or reverse
	
			// time table
			drcInflection << "time" << i_loopNumber << " <- c(";
			for (int i = 0; i < (i_bamFiles - 1); i++)  	
				drcInflection << turnoverTimes[i] << ",";	
			drcInflection << turnoverTimes[(i_bamFiles - 1)] << ")\n";

			// depth table
			drcInflection << "depth" << i_loopNumber << " <- c(";
			for (int i = 0; i < (i_bamFiles - 1); i++)  	
				drcInflection << turnoverArr[i] << ",";	
			drcInflection << turnoverArr[(i_bamFiles - 1)] << ")\n";

			drcInflection << drcR(i_loopNumber, std::string(""), b_L4flag);

		} else if ( ((b_typeArr == b_peak1Arr) || (b_typeArr == b_peak2Arr)) && (b_model == true) ) {	//peak1 or peak2

			// forward		
			int i_forwardEnd = i_maxIndex;
			for (int i = i_maxIndex; i < i_bamFiles; i++) {
				if (turnoverArr[i_maxIndex]* d_satthresh < turnoverArr[i]) {
					i_forwardEnd = i;
				}
			}

			// time table
			drcInflection << "time" << i_loopNumber << ".1 <- c(";
			for (int i = 0; i < i_forwardEnd; i++)  	
				drcInflection << turnoverTimes[i] << ",";	
			drcInflection << turnoverTimes[i_forwardEnd] << ")\n";

			// depth table
			drcInflection << "depth" << i_loopNumber << ".1 <- c(";
			for (int i = 0; i < i_forwardEnd; i++)  	
				drcInflection << turnoverArr[i] << ",";	
			drcInflection << turnoverArr[i_forwardEnd] << ")\n";

			drcInflection << drcR(i_loopNumber, std::string("1"), b_L4flag);

			// reverse		
			int i_reverseStart = i_maxIndex;
			for (int i = i_maxIndex; i > 0; i--) {
				if (turnoverArr[i_maxIndex]*d_satthresh < turnoverArr[i]) {
					i_reverseStart = i;
				}
			}

			// time table
			drcInflection << "time" << i_loopNumber << ".2 <- c(";
			for (int i = i_reverseStart; i < i_bamFiles-1; i++)  	
				drcInflection << turnoverTimes[i] << ",";	
			drcInflection << turnoverTimes[i_bamFiles-1] << ")\n";

			// depth table
			drcInflection << "depth" << i_loopNumber << ".2 <- c(";
			for (int i = i_reverseStart; i < i_bamFiles-1; i++)  	
				drcInflection << turnoverArr[i] << ",";	
			drcInflection << turnoverArr[i_bamFiles-1] << ")\n";

			drcInflection << drcR(i_loopNumber, std::string("2"), b_L4flag);


		} else if ( ((b_typeArr == b_dip1Arr) || (b_typeArr == b_dip2Arr)) && (b_model == true) ) { 	//dip1 or dip2

			// reverse		
			int i_reverseEnd = i_minIndex;
			for (int i = i_minIndex; i < i_bamFiles; i++) {
				if (turnoverArr[i]* d_satthresh < turnoverArr[i_minIndex]) {
					i_reverseEnd = i;
				}
			}

			// time table
			drcInflection << "time" << i_loopNumber << ".1 <- c(";
			for (int i = 0; i < i_reverseEnd; i++)  	
				drcInflection << turnoverTimes[i] << ",";	
			drcInflection << turnoverTimes[i_reverseEnd] << ")\n";

			// depth table
			drcInflection << "depth" << i_loopNumber << ".1 <- c(";
			for (int i = 0; i < i_reverseEnd; i++)  	
				drcInflection << turnoverArr[i] << ",";	
			drcInflection << turnoverArr[i_reverseEnd] << ")\n";

			drcInflection << drcR(i_loopNumber, std::string("1"), b_L4flag);

			// forward		
			int i_forwardStart = i_minIndex;
			for (int i = i_minIndex; i > 0; i--) {
				if (turnoverArr[i]* d_satthresh < turnoverArr[i_maxIndex]) {
					i_forwardStart = i;
				}
			}

			// time table
			drcInflection << "time" << i_loopNumber << ".2 <- c(";
			for (int i = i_forwardStart; i < i_bamFiles-1; i++)  	
				drcInflection << turnoverTimes[i] << ",";	
			drcInflection << turnoverTimes[i_bamFiles-1] << ")\n";

			// depth table
			drcInflection << "depth" << i_loopNumber << ".2 <- c(";
			for (int i = i_forwardStart; i < i_bamFiles-1; i++)  	
				drcInflection << turnoverArr[i] << ",";	
			drcInflection << turnoverArr[i_bamFiles-1] << ")\n";

			drcInflection << drcR(i_loopNumber, std::string("2"), b_L4flag);

		} else {											//undefined - model anyway
			// time table
			drcInflection << "time" << i_loopNumber << " <- c(";
			for (int i = 0; i < (i_bamFiles - 1); i++)  	
				drcInflection << turnoverTimes[i] << ",";	
			drcInflection << turnoverTimes[(i_bamFiles - 1)] << ")\n";

			// depth table
			drcInflection << "depth" << i_loopNumber << " <- c(";
			for (int i = 0; i < (i_bamFiles - 1); i++)  	
				drcInflection << turnoverArr[i] << ",";	
			drcInflection << turnoverArr[(i_bamFiles - 1)] << ")\n";

			drcInflection << drcR(i_loopNumber, std::string(""), b_L4flag);
		}
	}

}
