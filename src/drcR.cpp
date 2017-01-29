#include <fstream>
#include <iostream>
#include <string>
#include <regex>
#include <vector>
using namespace std;

std::string drcR(int i_loopNumber, std::string s_peakExtention)  
{

	// IF STATEMENT
	/** SUMMARY:
	let:
		L.5      = regular 5 parameter logistic (L5) 
		FL.5     = L5 with lower asyptote forced to 0
		FLMIN.5  = L5 with lower asyptote forced to minimum value in depth vector
		rL.5     = same as L5 except upper asymptote fixed to maximum value in depth vector
		rFL.5    = same as FL.5 except upper asymptote fixed to maximum value in depth vector
		rFLMIN.5 = same as FLMIN.5 except upper asymptote fixed to maximum value in depth vector

	if ( (L.5 lower asyptote < 0) && (L.5 = forward) && (FL.5 inflection point < (i_Mult*minTime) ) { FL.5     }
	else if ( (L.5 lower asyptote >= 0) && (L.5 = forward) ) 							{ L.5      }
	else if ( (L.5 lower asyptote < 0) && (L.5 = forward) ) 							{ FLMIN.5  }
	else if ( (rL.5 lower asyptote < 0) && (rFL.5 inflection point < (i_Mult*minTime) )             { rFL.5    }
	else if (rL.5 lower asyptote >= 0)                                                              { rL.5     }
	else 				                                                                        { rFLMIN.5 } 
	**/

	std::string s_rCommand = ""; // grow this and return
	std::string s_LN = to_string(i_loopNumber); 
	std::string s_EX; // s_EX = sextention 
	if (s_peakExtention == "") { s_EX = s_LN; }
	else { s_EX = s_LN + "." + s_peakExtention; }
	int i_Mult = 5; // multiplier for inflection point cut off limit 
	std::string s_Mult = to_string(i_Mult); 

	s_rCommand += "peak" + s_EX + " <- data.frame(time" + s_EX + ",depth" + s_EX + ")\n";
	s_rCommand += "print(paste(\"peak" + s_EX + "\"));\n";
	s_rCommand += "suppressMessages(library(drc))\n";

	// L4 function
	s_rCommand += "try(peak" + s_EX + ".L.5 <- drm(depth" + s_EX + " ~ time" + s_EX + ", data = peak" + s_EX + ", fct = L.5()))\n";
	s_rCommand += "try(peak" + s_EX + ".FL.5 <- drm(depth" + s_EX + " ~ time" + s_EX + ", data = peak" + s_EX + ", fct = L.5(fixed=c(NA,0,NA,NA,NA))))\n";
	s_rCommand += "try(peak" + s_EX + ".FLMIN.5 <- drm(depth" + s_EX + " ~ time" + s_EX + ", data = peak" + s_EX + ", fct = L.5(fixed=c(NA,min(depth" + s_EX + "),NA,NA,NA))))\n";
	s_rCommand += "try(Rpeak" + s_EX + ".L.5 <- drm(depth" + s_EX + " ~ time" + s_EX + ", data = peak" + s_EX + ", fct = L.5(fixed=c(NA,NA,max(depth" + s_EX + "),NA,NA))))\n"; 
	s_rCommand += "try(Rpeak" + s_EX + ".FL.5 <- drm(depth" + s_EX + " ~ time" + s_EX + ", data = peak" + s_EX + ", fct = L.5(fixed=c(NA,0,max(depth" + s_EX + "),NA,NA))))\n"; 
	s_rCommand += "try(Rpeak" + s_EX + ".FLMIN.5 <- drm(depth" + s_EX + " ~ time" + s_EX + ", data = peak" + s_EX + ", fct = L.5(fixed=c(NA,min(depth" + s_EX + "),max(depth" + s_EX + "),NA,NA))))\n"; 

	// IF STATEMENT

	s_rCommand += "if ( (exists(\"peak" + s_EX + ".L.5\")) & (exists(\"peak" + s_EX + ".FL.5\")) & (exists(\"peak" + s_EX + ".FLMIN.5\")) & (exists(\"Rpeak" + s_EX +
		".L.5\")) & (exists(\"Rpeak" + s_EX + ".FL.5\")) & (exists(\"Rpeak" + s_EX + ".FLMIN.5\")) ) {\n";
	// FL.5
	s_rCommand += "  if ( (peak" + s_EX + ".L.5$parmMat[2] < 0) & (peak" + s_EX + ".L.5$parmMat[1] < 0) & (peak" + s_EX + ".FL.5$parmMat[3] < (" + s_Mult + "*max(time" + s_EX + "))) ) {\n";
	s_rCommand += "    tryCatch(suppressWarnings(summary(peak" + s_EX + ".FL.5)), print(paste(\"c:(Intercept) FORCED\")), error = function(e) {print(paste(\"e:(Intercept) ERROR\")); " +
		"print(paste(\"d:(Intercept) ERROR\")); print(paste(\"c:(Intercept) ERROR\")); print(paste(\"b:(Intercept) ERROR\")); print(paste(\"f:(Intercept) ERROR\")); " +
		"print(paste(\" NaN (-1 degrees of freedom)\")); NaN})\n";
	// L.5
	s_rCommand += "  } else if ( (peak" + s_EX + ".L.5$parmMat[2] >= 0) & (peak" + s_EX + ".L.5$parmMat[1] < 0) ) {\n";
	s_rCommand += "    tryCatch(suppressWarnings(summary(peak" + s_EX + ".L.5)),error = function(e) {print(paste(\"e:(Intercept) ERROR\")); print(paste(\"d:(Intercept) ERROR\")); " + 
		"print(paste(\"c:(Intercept) ERROR\")); print(paste(\"b:(Intercept) ERROR\")); print(paste(\"f:(Intercept) ERROR\")); print(paste(\" NaN (-1 degrees of freedom)\")); NaN})\n";
	// FLMIN.5
	s_rCommand += "  } else if ( (peak" + s_EX + ".L.5$parmMat[2] < 0) & (peak" + s_EX + ".L.5$parmMat[1] < 0) ) {\n";
	s_rCommand += "    tryCatch(suppressWarnings(summary(peak" + s_EX + ".FLMIN.5)), cat(\"c:(Intercept)\", min(depth" + s_EX + "), \"\\n\"), error = function(e) {" +
		"print(paste(\"e:(Intercept) ERROR\")); print(paste(\"d:(Intercept) ERROR\")); print(paste(\"c:(Intercept) ERROR\")); print(paste(\"b:(Intercept) ERROR\")); " +
		"print(paste(\"f:(Intercept) ERROR\")); print(paste(\" NaN (-1 degrees of freedom)\")); NaN})\n";
	// rFL.5
	s_rCommand += "  } else if ( (Rpeak" + s_EX + ".L.5$parmMat[2] < 0) & (Rpeak" + s_EX + ".FL.5$parmMat[2] < (" + s_Mult + "*max(time" + s_EX + "))) ) {\n";
	s_rCommand += "    tryCatch(suppressWarnings(summary(Rpeak" + s_EX + ".FL.5)), print(paste(\"c:(Intercept) FORCED\")), cat(\"d:(Intercept)\", max(depth" + s_EX + "), \"\\n\"), " +
		"error = function(e) {print(paste(\"e:(Intercept) ERROR\")); print(paste(\"d:(Intercept) ERROR\")); print(paste(\"c:(Intercept) ERROR\")); " + 
		"print(paste(\"b:(Intercept) ERROR\")); print(paste(\"f:(Intercept) ERROR\")); print(paste(\" NaN (-1 degrees of freedom)\")); NaN})\n";
	// rL.5
	s_rCommand += "  } else if (Rpeak" + s_EX + ".L.5$parmMat[2] >= 0) {\n";
	s_rCommand += "    tryCatch(suppressWarnings(summary(Rpeak" + s_EX + ".L.5)), cat(\"d:(Intercept)\", max(depth" + s_EX + "), \"\\n\"), error = function(e) " +
		"{print(paste(\"e:(Intercept) ERROR\")); print(paste(\"d:(Intercept) ERROR\")); print(paste(\"c:(Intercept) ERROR\")); " + 
		"print(paste(\"b:(Intercept) ERROR\")); print(paste(\"f:(Intercept) ERROR\")); print(paste(\" NaN (-1 degrees of freedom)\")); NaN})\n";
	// rFLMIN.5
 	s_rCommand += "  } else {\n";
	s_rCommand += "    tryCatch(suppressWarnings(summary(Rpeak" + s_EX + ".FLMIN.5)), cat(\"d:(Intercept)\", max(depth" + s_EX + "), \"\\n\"), cat(\"c:(Intercept)\", min(depth" + 
		s_EX + "), \"\\n\"), error = function(e) {print(paste(\"e:(Intercept) ERROR\")); print(paste(\"d:(Intercept) ERROR\")); print(paste(\"c:(Intercept) ERROR\")); " +
		"print(paste(\"b:(Intercept) ERROR\")); print(paste(\"f:(Intercept) ERROR\")); print(paste(\" NaN (-1 degrees of freedom)\")); NaN})\n";
 	s_rCommand += "  }\n";

 	s_rCommand += "} else {\n";
 	s_rCommand += "  print(paste(\"e:(Intercept) ERROR\")); print(paste(\"d:(Intercept) ERROR\")); print(paste(\"c:(Intercept) ERROR\")); print(paste(\"b:(Intercept) ERROR\")); print(paste(\"f:(Intercept) ERROR\"))\n";
 	s_rCommand += "}\n";




	return s_rCommand;

}

