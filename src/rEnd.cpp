#include <iostream>
#include <fstream>
#include <string>
using namespace std;

void rEnd(std::string s_rPltsName,std::string s_rPltsPDF, bool b_genomeFlagCall, int i_bamRep, bool b_model, int i_increase, int i_decrease, bool b_linFlag) 
{
	std::ofstream rPlts;
	rPlts.open(s_rPltsName, std::ios::app);

	rPlts << "pdf(\"" << s_rPltsPDF << "\")\n";
	rPlts << "piebarLayout <- matrix(c(1,2), nrow = 2, byrow = TRUE)\n";
	rPlts << "suppressWarnings(multiplot(pie,bar, layout = piebarLayout))\n";
	rPlts << "suppressWarnings(print(heatmap))\n";
	// Average profiles
	rPlts << "multiplot(riseProfile,fallProfile, layout = piebarLayout)\n";

	if (b_model) { 
		rPlts << "modelLayout <- matrix(c(1,2,3,4,5,6), nrow = 3, byrow = TRUE)\n";
		rPlts << "multiplot(UpperA, LowerA, undefRiseProfile, undefFallProfile, hillProfile, valleyProfile, layout = modelLayout)\n";
	}

	if (i_bamRep > 1) {
		rPlts << "layout1 <- matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE)\n";
		rPlts << "suppressWarnings(multiplot(s1,s2,stdevBox, layout = layout1))\n";
	} 

	if (b_model) {
		rPlts << "scatLayout <- matrix(c(1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 3, 2, 2, 2, 4), nrow = 4, byrow = TRUE)\n";
		if ( (i_increase > 0) && (i_decrease > 0) ) {
			rPlts << "suppressWarnings(multiplot(forwardScat,forwardStackBarInf,forwardStackBarHC,legFor, layout = scatLayout))\n";
			rPlts << "suppressWarnings(multiplot(reverseScat,reverseStackBarInf,reverseStackBarHC,legRev, layout = scatLayout))\n";
		} else if (i_increase > 0) {
			rPlts << "multiplot(forwardScat,forwardStackBarInf,forwardStackBarHC,legFor, layout = scatLayout)\n";
		} else if (i_decrease > 0) {
			rPlts << "multiplot(reverseScat,reverseStackBarInf,reverseStackBarHC,legRev, layout = scatLayout)\n";
		} else {}
	} else {
		rPlts << "scatLayout <- matrix(c(1, 1, 1, 3, 1, 1, 1, 3, 1, 1, 1, 3, 2, 2, 2, 0), nrow = 4, byrow = TRUE)\n";
		if ( (i_increase > 0) && (i_decrease > 0) ) {
			rPlts << "suppressWarnings(multiplot(forwardScat,forwardStackBarInf,forwardStackBarHC, layout = scatLayout))\n";
			rPlts << "suppressWarnings(multiplot(reverseScat,reverseStackBarInf,reverseStackBarHC, layout = scatLayout))\n";
		} else if (i_increase > 0) {
			rPlts << "multiplot(forwardScat,forwardStackBarInf,forwardStackBarHC, layout = scatLayout)\n";
		} else if (i_decrease > 0) {
				rPlts << "multiplot(reverseScat,reverseStackBarInf,reverseStackBarHC, layout = scatLayout)\n";
		} else {}
	}

	if (b_linFlag) {
		rPlts << "linRegLayout <- matrix(c(1,2,3,4), nrow = 2, byrow = TRUE)\n";
		rPlts << "suppressWarnings(multiplot(dCov,sigRes,slope,linRes, layout = linRegLayout))\n";
	} else {
		rPlts << "suppressWarnings(multiplot(dCov,sigRes, layout = piebarLayout))\n";
	}

	if ( b_genomeFlagCall == true ) {
		if ( (i_increase > 0) && (i_decrease > 0) ) {
			rPlts << "boxLayout <- matrix(c(0,1,1,1,1,0,0,2,2,2,2,0), nrow = 2, byrow = TRUE)\n";
			rPlts << "suppressWarnings(multiplot(boxForward,boxReverse, layout = boxLayout))\n";
			rPlts << "suppressWarnings(print(ideForward))\n";
			rPlts << "suppressWarnings(print(ideReverse))\n";
		} else if (i_increase > 0) {
			rPlts << "suppressWarnings(print(boxForward))\n";
			rPlts << "suppressWarnings(print(ideForward))\n";
		} else if (i_decrease > 0) {
			rPlts << "suppressWarnings(print(boxReverse))\n";
			rPlts << "suppressWarnings(print(ideReverse))\n";
		} else {}
	}

	// OTHER CHARTS
	rPlts << "invisible(dev.off())\n";
	rPlts.close();

	// run rthon script
	std::string s_rscript = "Rscript " + s_rPltsName; // temp file
	
	try {
		std::system(s_rscript.c_str());
	}
	catch(...){
		std::cerr << "Cannot access R" << std::endl;
		std::exit(0);
	}
}
