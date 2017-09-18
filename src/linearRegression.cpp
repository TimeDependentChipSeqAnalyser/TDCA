#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

std::string exec(const char* cmd);

void linearRegression(int i_bamFiles, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, int turnoverTimes[], std::string s_name) {

	std::string s_rLinRegName = s_name + ".tdca.LinReg.R";			// Temporary file
	std::string s_regexLinRegName = s_name + ".tdca.LinReg.regex.txt";	// Temporary file

	std::string s_yint = "yint.tdca.txt";						// Temporary file
	std::string s_residuals = "linRegRes.tdca.txt";					// Temporary file

	std::string s_txtLinRegName = s_name + ".tdca.LinReg.txt"; 			// Output file

	// R script
	std::ofstream rLinReg;
	rLinReg.open (s_rLinRegName);
	rLinReg << "time = c(";
	for (int i = 0; i < i_bamFiles-1; i++) {
		rLinReg << turnoverTimes[i] << ",";
	}
	rLinReg << turnoverTimes[i_bamFiles-1] << ")\n";

	//dataArray[row][col]
	for (int i = 1; i < i_peakNumber+1; i++) {
		rLinReg << "coverage" << i << " = c(";

		for (int j = 3; j < i_bamFiles+2; j++) {
			rLinReg << dataArray[i][j] << ",";
		}
		rLinReg << dataArray[i][i_bamFiles+2] << ")\n";

		rLinReg << "peak" << i << " <- lm(coverage" << i << " ~ time)\n";
		rLinReg << "summary(peak" << i << ")\n";
	}

	rLinReg.close();

	// run r script
	std::string s_rScript = "Rscript " + s_rLinRegName + " > " + s_regexLinRegName; 
	try {
		std::system(s_rScript.c_str());
	}
	catch(...){
		std::cerr << "Cannot access R" << std::endl;
		std::exit(0);
	}

	// REGEX
	// y-int
	std::string s_regex = "grep '(Intercept)' " + s_regexLinRegName + " | grep -o '[^[:space:]]*[[:space:]]*[^[:space:]]*' |" +
		"grep '(Intercept)' | sed -n 's/(Intercept)//p' | grep -o '[^[:space:]]*' > yint.tdca.txt";
	try {
		std::system(s_regex.c_str());
	} catch(...){
		std::cerr << "Linear regression y-int script failed.";
		std::exit(0);
	}

	// residuals
	s_regex = "sed -n 's/Residual standard error: //p' " + s_regexLinRegName + " | sed 's/on.*//' > linRegRes.tdca.txt";
	try {
		std::system(s_regex.c_str());
	} catch(...){
		std::cerr << "Linear regression residuals script failed.";
		std::exit(0);
	}

	// Write y-int and residuals to vectors
	double d_yint[i_peakNumber];
	double d_residuals[i_peakNumber];
	double d_slope[i_peakNumber];

	std::string s_line;

	std::ifstream yintFile ("yint.tdca.txt");
	if (yintFile.is_open()) {
		int i_count = 0;
		while ( std::getline (yintFile, s_line) ) {
			d_yint[i_count] = stod(s_line);
			i_count++;
		}
		yintFile.close();
	} else { 
		std::cout << "Unable to read linear regression y-int file. Program terminated."; 
		std::exit(0);
	}

	std::ifstream resFile ("linRegRes.tdca.txt");
	if (resFile.is_open()) {
		int i_count = 0;
		while ( std::getline (resFile, s_line) ) {
			d_residuals[i_count] = stod(s_line);
			i_count++;
		}
		resFile.close();
	} else { 
		std::cout << "Unable to read linear regression residuals file. Program terminated."; 
		std::exit(0);
	}


	for (int i = 0; i < i_peakNumber; i++) {
		// use coverage at second time point and value of second time point for 'y' and 'x' respectively
		d_slope[i] = (stod(dataArray[i+1][4])- d_yint[i]) / double(turnoverTimes[1]);
	}

	// PRINT TO OUTPUT FILE
	std::ofstream linRegOut;
	linRegOut.open (s_txtLinRegName);

	//dataArray[row][col]
	// Header first:
	for (int j = 0; j < i_bamFiles+3; j++) {
		linRegOut << dataArray[0][j] << "\t";
	}
	linRegOut << "Y intercept\tSlope\tResiduals\n";


	// Contents:
	for (int i = 1; i < i_peakNumber; i++) {
		for (int j = 0; j < i_bamFiles+3; j++) {
			linRegOut << dataArray[i][j] << "\t";
		}
		linRegOut << d_yint[i-1] << "\t" << d_slope[i-1] << "\t" << d_residuals[i-1] << "\n";
	}

	for (int j = 0; j < i_bamFiles+3; j++) {
		linRegOut << dataArray[i_peakNumber][j] << "\t";
	}

	linRegOut << d_yint[i_peakNumber-1] << "\t" << d_slope[i_peakNumber-1] << "\t" << d_residuals[i_peakNumber-1];

	linRegOut.close();




	// Remove temporary files
	std::remove(std::string(s_rLinRegName).c_str());
	std::remove(std::string(s_regexLinRegName).c_str());
	std::remove(std::string(s_yint).c_str());
	std::remove(std::string(s_residuals).c_str());
}
