#include <fstream>
#include <iostream>
#include <string>
#include <regex>
#include <vector>
using namespace std;

std::string drcRpoisson(int i_loopNumber, std::string s_peakExtention)  
{

	std::string s_rCommand = ""; // grow this and return
	std::string s_LN = to_string(i_loopNumber); 
	std::string s_EX; // s_EX = sextention 
	if (s_peakExtention == "") { s_EX = s_LN; }
	else { s_EX = s_LN + "." + s_peakExtention; }

	s_rCommand += "peak" + s_EX + " <- data.frame(time" + s_EX + ",depth" + s_EX + ")\n";
	s_rCommand += "print(paste(\"peak" + s_EX + "\"));\n";
	s_rCommand += "suppressMessages(library(drc))\n";


	// poisson L3 function
	s_rCommand += "try(peak" + s_EX + ".L.3 <- drm(depth" + s_EX + " ~ time" + s_EX + ", data = peak" + s_EX + ", fct = L.3(), type = \"Poisson\"))\n";


	// IF STATEMENT
	s_rCommand += "if (exists(\"peak" + s_EX + ".L.3\") ) {\n"; // if poisson model exists
	s_rCommand += "    tryCatch(suppressWarnings(summary(peak" + s_EX + ".L.3)), print(paste(\"c:(Intercept) 0\")), print(paste(\"f:(Intercept) 1\")), print(paste(\" NaN (-1 degrees of freedom)\")), " + 
		"error = function(e) {print(paste(\"e:(Intercept) ERROR\")); print(paste(\"d:(Intercept) ERROR\")); print(paste(\"c:(Intercept) ERROR\")); " + 
		"print(paste(\"b:(Intercept) ERROR\")); print(paste(\"f:(Intercept) ERROR\")); print(paste(\" NaN (-1 degrees of freedom)\")); NaN})\n";

	s_rCommand += "} else {\n";
	s_rCommand += "  print(paste(\"e:(Intercept) ERROR\")); print(paste(\"d:(Intercept) ERROR\")); print(paste(\"c:(Intercept) ERROR\")); print(paste(\"b:(Intercept) ERROR\")); print(paste(\"f:(Intercept) ERROR\")); print(paste(\" NaN (-1 degrees of freedom)\"))\n";
	s_rCommand += "}\n";


	return s_rCommand;

}
