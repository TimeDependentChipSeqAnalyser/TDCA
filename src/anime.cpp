#include <iostream>
#include <fstream>
#include <string>
#include <vector>
//#include <regex>
#include <algorithm>
#include <cmath>
#include <math.h> 
using namespace std;

std::string exec(const char* cmd);
void r3Dadjusted(std::string s_genelistFile, std::string s_rGenePltsName, int i_bamFiles, int turnoverTimes[], int i_validgenes, std::vector<double> &normGeneVector, std::vector<std::string> &validGeneVector);
std::string strand(std::string s_gene, std::string s_genelistFile);
int exonNumber(std::string s_geneName, std::string s_genelistFile);
std::string exonStartString(std::string s_geneName, std::string s_genelistFile, int i_exons);
std::string exonEndString(std::string s_geneName, std::string s_genelistFile, int i_exons);
int exons(std::string s_exonStartString);
std::string exonStringSnipper(std::string s_exonStartString);
int exons(std::string s_exonEndString);
std::string exonStringSnipper(std::string s_exonEndString);
double startEnd(std::string s_geneName, std::string s_genelistFile, std::string);

void anime(int i_bamFiles, int i_validgenes, std::vector<double> &normGeneVector, int turnoverTimes[], std::string s_name, std::vector<std::string> &validGeneVector, std::string s_genelistFile) {

	double d_barCount = 1200; // increase for better resolution

	// d_smallestTimeGap and d_timeGaps[] are constant for all genes
	double d_smallestTimeGap = turnoverTimes[1] - turnoverTimes[0];
	for (int i = 1; i < i_bamFiles-1; i++) {
		int i_temp = turnoverTimes[i+1] - turnoverTimes[i];
		if (i_temp < d_smallestTimeGap) { d_smallestTimeGap = i_temp; }
	}
	double d_timeGaps[i_bamFiles-1];
	for (int i = 0; i < i_bamFiles-1; i++) {
		d_timeGaps[i] = d_smallestTimeGap / (turnoverTimes[i+1] - turnoverTimes[i]);
	}


	// loop through all genes
	for (int i = 0; i < i_validgenes; i++) {

		ofstream jScript;
		std::string s_fileName = s_name + "." + validGeneVector[i] + ".html";
		jScript.open (s_fileName);

		jScript << "<!DOCTYPE html>\n";
		jScript << "<html>\n";
		jScript << "<style>\n";
		jScript << "#myContainer { width: 1350px; height: 625px; position: absolute; background: white; }\n";

		// draw exons
		std::string s_geneName = validGeneVector[i];
		int i_exons = exonNumber(s_geneName, s_genelistFile);
		std::string s_exonStartString = exonStartString(s_geneName, s_genelistFile, i_exons);
		std::string s_exonEndString = exonEndString(s_geneName, s_genelistFile, i_exons);

		std::vector<double> exonStartVector(0);
		std::vector<double> exonEndVector(0);

		// start exon vector
		for (int k = 0; k < i_exons-1; k++) { 
			int i_anExon = exons(s_exonStartString);
			exonStartVector.push_back(i_anExon);
			s_exonStartString = exonStringSnipper(s_exonStartString);
		}	
		exonStartVector.push_back(stoi(s_exonStartString));

		// end exon vector
		for (int k = 0; k < i_exons-1; k++) { 
			int i_anExon = exons(s_exonEndString);
			exonEndVector.push_back(i_anExon);
			s_exonEndString = exonStringSnipper(s_exonEndString);
		}	
		exonEndVector.push_back(stoi(s_exonEndString));

		std::string s_strand = strand(s_geneName, s_genelistFile);

		double d_start = startEnd(s_geneName, s_genelistFile, std::string ("start"));
		double d_end = startEnd(s_geneName, s_genelistFile, std::string ("end"));

		double d_length = d_end-d_start;
		if (s_strand == "+") { // strand is positive
			for (int k = 0; k < i_exons; k++) { 
				double d_boxStart = (exonStartVector[k] - d_start) / d_length; 
				double d_boxEnd = (exonEndVector[k] - d_start) / d_length;
				jScript << "#exon" << k+1 << " { width: " << (d_boxEnd-d_boxStart)*800 << "px; height: 20px; position: absolute; left: " << 300+d_boxStart*800 << "px; top: 550px; background-color: black; }\n";
	
			} 
		} else { // strand is negtive
			for (int k = i_exons-1; k > -1; k--) { 
				double d_boxStart = (d_end - exonEndVector[k]) / d_length; 
				double d_boxEnd = (d_end - exonStartVector[k]) / d_length;
				jScript << "#exon" << k+1 << " { width: " << (d_boxEnd-d_boxStart)*800 << "px; height: 20px; position: absolute; left: " << 300+d_boxStart*800 << "px; top: 550px; background-color: black; }\n";	
			}
		}
		// ===================================================== //
		// XXX HTML BABY XXX HTML BABY XXX SOOOOOO GOOOOOOOD XXX //
		// ===================================================== //
		jScript << "#gene { width: 800px; height: 2px; position: absolute; left: 300px; top: 560px; background-color: black; }\n";
		jScript << "#name { width: 200px; top: 10px; height: 50px; left: 650px; position: absolute; background: white; font-size: 18px; text-align: center; }\n";
		jScript << "#tb { width: 200px; top: 0px; height: 50px; left: 120px; position: absolute; background: white; }\n";
		jScript << "#hl { width: 1202px; height: 2px; position: absolute; left: 100px; top: 525px; background-color: black; }\n";
		jScript << "#ht1 { width: 4px; height: 2px; position: absolute; left: 96px; top: 525px; background-color: black; }\n";
		jScript << "#ht2 { width: 4px; height: 2px; position: absolute; left: 96px; top: 400px; background-color: black; }\n";
		jScript << "#ht3 { width: 4px; height: 2px; position: absolute; left: 96px; top: 275px; background-color: black; }\n";
		jScript << "#ht4 { width: 4px; height: 2px; position: absolute; left: 96px; top: 125px; background-color: black; } \n";
		jScript << "#ht5 { width: 4px; height: 2px; position: absolute; left: 96px; top: 25px; background-color: black; } \n";
		jScript << "#depth1 { width: 50px; height: 2px; position: absolute; left: 50px; top: 517px; font-size: 13.5px; text-align: left; } \n";
		jScript << "#depth2 { width: 50px; height: 2px; position: absolute; left: 50px; top: 392px; font-size: 13.5px; text-align: left; } \n";
		jScript << "#depth3 { width: 50px; height: 2px; position: absolute; left: 50px; top: 267px; font-size: 13.5px; text-align: left; } \n";
		jScript << "#depth4 { width: 50px; height: 2px; position: absolute; left: 50px; top: 117px; font-size: 13.5px; text-align: left; } \n";
		jScript << "#depth5 { width: 50px; height: 2px; position: absolute; left: 50px; top: 17px; font-size: 13.5px; text-align: left; } \n";
		jScript << "#vl { width: 2px; height: 500px; position: absolute; left: 100px; top: 25px; background-color: black; }\n";
		jScript << "div.absolute1 { position: absolute; top: 0px; left: 120px; width: 50px; height: 50px; }\n";
		jScript << "div.absolute2 { position: absolute; top: 0px; left: 180px; width: 200px; height: 50px; }\n";
		jScript << "#timebar { width: 5px; height: 10px; position: absolute; left: 102px; top: 580px; background-color: red; }\n";
		jScript << "#timebox { width: 1200px; height: 10px; position: absolute; left: 100px; top: 578px; border-style: solid; border-width: 2px; border-color: black; }\n";
		jScript << "#tp1 { width: 45px; height: 2px; position: absolute; left: 73px; top: 600px; transform: rotate(-90deg); font-size: 13.5px; text-align: left; }\n"; // first
		jScript << "#tp2 { width: 45px; height: 2px; position: absolute; left: 1273px; top: 600px; transform: rotate(-90deg); font-size: 13.5px; text-align: left; }\n";
		jScript << "#timetick1 { width: 2px; height: 4px; position: absolute; left: 100px; top: 592px; background-color: black; }\n"; // last
		jScript << "#timetick2 { width: 2px; height: 4px; position: absolute; left: 1302px; top: 592px; background-color: black; }\n";
		for (int j = 1; j < i_bamFiles-1; j++) {
			double d_d = (double(turnoverTimes[j])/double(turnoverTimes[i_bamFiles-1]))*1200;
			jScript << "#tp" << j+2 << " { width: 45px; height: 2px; position: absolute; left: " << d_d + 73 << "px; top: 600px; transform: rotate(-90deg); font-size: 13.5px; text-align: left; }\n";
			jScript << "#timetick" << j+2 << " { width: 2px; height: 4px; position: absolute; left: " << d_d + 102 << "px; top: 592px; background-color: black; }\n"; // last
		}

		std::vector<double> depthVector;
		// Print out normGeneVector contents
		for (int k = 0; k < i_bamFiles-1; k++) {
			for (int j = 0; j < d_barCount; j++) { 
				double d = normGeneVector[j+(k*d_barCount)+(i*i_bamFiles*d_barCount)];
				depthVector.push_back(d);
			}
		}
		for (int j = 0; j < d_barCount-1; j++) { 
			double d = normGeneVector[j+((i_bamFiles-1)*d_barCount)+(i*i_bamFiles*d_barCount)];  
			depthVector.push_back(d);
		}
		double d = normGeneVector[d_barCount-1+((i_bamFiles-1)*d_barCount)+(i*i_bamFiles*d_barCount)];
		depthVector.push_back(d);

		// figure out y axis max/min
		double d_max = 0;
		double d_min = depthVector[0];
		for (int k = 0; k < depthVector.size(); k++) {			
			double d_temp = depthVector[k];
			if (d_temp > d_max) { d_max = d_temp; }
			if (d_temp < d_min) { d_min = d_temp; }
		}

		// figure out max y axis range
		double d_range = 0;
		for (int j = 0; j < d_barCount; j++) { 
			std::vector<double> d_pointVector;
			for (int k = 0; k < i_bamFiles; k++) {
				double d = normGeneVector[j+(k*d_barCount)+(i*i_bamFiles*d_barCount)];
				d_pointVector.push_back(d);
			}
			double d_pointMax = 0;
			double d_pointMin = d_pointVector[0];
			for (int k = 0; k < d_pointVector.size(); k++) {
				double d_temp = d_pointVector[k];
				if (d_temp > d_pointMax) { d_pointMax = d_temp; }
				if (d_temp < d_pointMin) { d_pointMin = d_temp; }
			}
			if ( (d_pointMax-d_pointMin) > d_range ) { d_range =  (d_pointMax-d_pointMin); }
			d_pointVector.clear();

		}

		// Print bars
		for (int j = 0; j < d_barCount; j++) {
			jScript << "#bar" << j+1 << " { width: 1px; height: " << (depthVector[j]-d_min)*500/d_range << "px; position: absolute; left: " <<
				102+j*1 << "px; top: " << 525 - ((depthVector[j]-d_min)*500/d_range) << "px; background-color: grey; }\n";
		} 

		jScript << "</style>\n";
		jScript << "<body>\n";
		jScript << "<p>\n";
		jScript << "<button onclick=\"myMove()\">Start</button> <!Change \"Start\" to gene name>\n";
		jScript << "<button onclick=\"reset()\">Reset</button> <!Change \"Start\" to gene name>\n";
		jScript << "</p>\n";
		jScript << "<div id =\"tb\"></div>\n";
		jScript << "<div id =\"myContainer\">\n";
		jScript << "	<div id =\"gene\"></div>\n";
		for (int j = 0; j < i_exons; j++) { 
			jScript << "	<div id =\"exon" << j+1 << "\"></div>\n";
		}
		jScript << "	<div id =\"hl\"></div>\n";
		jScript << "	<div id =\"ht1\"></div>\n";
		jScript << "	<div id =\"ht2\"></div>\n";
		jScript << "	<div id =\"ht3\"></div>\n";
		jScript << "	<div id =\"ht4\"></div>\n";
		jScript << "	<div id =\"ht5\"></div>\n";
		jScript << "	<div id =\"vl\"></div>\n";
		jScript << "	<div id =\"depth1\"><b>" << d_min << "</b></div>\n";
		jScript << "	<div id =\"depth2\"><b>" << (d_min + 0.25*d_range) << "</b></div>\n";
		jScript << "	<div id =\"depth3\"><b>" << d_min + 0.50*d_range << "</b></div>\n";
		jScript << "	<div id =\"depth4\"><b>" << d_min + 0.75*d_range << "</b></div>\n";
		jScript << "	<div id =\"depth5\"><b>" << d_max << "</b></div>\n";
		jScript << "	<div id =\"timebar\"></div>\n";
		jScript << "	<div id =\"timebox\"></div>\n";
		jScript << "	<div id =\"tp1\"><b>" << turnoverTimes[0] << "</b></div>\n";
		jScript << "	<div id =\"tp2\"><b>" << turnoverTimes[i_bamFiles-1] << "</b></div>\n";
		jScript << "	<div id =\"timetick1\"></div>\n";
		jScript << "	<div id =\"timetick2\"></div>\n";
		for (int j = 1; j < i_bamFiles-1; j++) {
			jScript << "	<div id =\"tp" << j+2 << "\"><b>" << turnoverTimes[j] << "</b></div>\n";
			jScript << "	<div id =\"timetick" << j+2 << "\"></div>\n";
		}

		for (int j = 0; j < d_barCount; j++) {
			jScript << "	<div id =\"bar" << j+1 << "\"></div>\n";
		} 
		jScript << "</div>\n";
		jScript << "<div class=\"absolute1\">\n";
		jScript << "	<p id=\"time\">TIME:</p>\n";
		jScript << "</div>\n";
		jScript << "<div class=\"absolute2\">\n";
		jScript << "	<p id=\"timeCount\">" << turnoverTimes[0] << "</p>\n";
		jScript << "</div>\n";
		jScript << "<div id =\"name\">" << validGeneVector[i] << "</div>\n";

		// ===================================================== //
		// XXX JAVASCIPT BAAAAAAAABYYY XXX SOOOOOO GOOOOOOOD XXX //
		// ===================================================== //

		jScript << "<script>\n";

		jScript << "var allBars = [";
		for (int j = 0; j < d_barCount-1; j++) {
			jScript << "\"bar" << j+1 << "\", ";
		}
		jScript << "\"bar" << d_barCount << "\"];\n";

		jScript << "var baseLine = [\n";
		for (int k = 0; k < i_bamFiles-1; k++) {
			jScript << "  [" << ( (depthVector[0+k*d_barCount]-d_min)*500/d_range ) + 1 << ", ";
			for (int j = 1; j < d_barCount-1; j++) {
				jScript << ( (depthVector[j+k*d_barCount]-d_min)*500/d_range ) + 1 << ", ";
			}
			jScript << ( (depthVector[(d_barCount-1)+k*d_barCount]-d_min)*500/d_range ) + 1 << "],\n";
		}
		jScript << "];\n";

		jScript << "var time = document.getElementById(\"timebar\");\n";

		jScript << "function myMove() {\n";
		jScript << "	var pos = [\n";
		for (int k = 0; k < i_bamFiles-1; k++) {
			jScript << "	  [0, ";
			for (int j = 1; j < d_barCount-1; j++) {
				jScript << "0, ";
			}
			jScript << "0],\n";
		}
		jScript << "	];\n";

		double d_1 = double(turnoverTimes[i_bamFiles-1]);
		double d_2 = double(turnoverTimes[1]);
		double d_3 = double(i_bamFiles);
		double d_adjuster = d_1 / (d_2*d_3);

		jScript << "	var posJumps = [\n";
		for (int k = 0; k < i_bamFiles-1; k++) {

			jScript << "	  [" << ( (depthVector[0+((k+1)*d_barCount)]-d_min) - (depthVector[0+((k)*d_barCount)]-d_min) ) / d_range * d_timeGaps[k] * d_adjuster  << ", "; 
			for (int j = 1; j < d_barCount-1; j++) {
				jScript << ( (depthVector[j+((k+1)*d_barCount)]-d_min) - (depthVector[j+((k)*d_barCount)]-d_min) ) / d_range * d_timeGaps[k] * d_adjuster << ", ";
			}			
			jScript << ( (depthVector[(d_barCount-1)+((k+1)*d_barCount)]-d_min) - (depthVector[(d_barCount-1)+((k)*d_barCount)]-d_min) ) / d_range * d_timeGaps[k] * d_adjuster << "],\n";

		}
		jScript << "];\n";

		jScript << "	var timePoints = " << i_bamFiles << ";\n";
		jScript << "	var id = setInterval(frame, 10); // speed\n";
		jScript << "	var counter = 0;\n";
		jScript << "	var speedup = 1;\n";
		jScript << "	var s_startTime = document.getElementById(\"timeCount\").innerHTML;\n";
		jScript << "	var i_startTime = parseInt(s_startTime);\n";
		jScript << "	var i_endTime = " << turnoverTimes[i_bamFiles-1] << "; // VARIABLE - relative units\n";

		jScript << "	function frame() {\n";
		jScript << "		for (i = 0; i < allBars.length; i++) {\n";
		jScript << "			time.style.left = 102 + (counter / (500*allBars.length*timePoints) * 1196 ) + 'px';\n";
		jScript << "			var elem = document.getElementById(allBars[i]);\n";
		jScript << "			if (counter == (500*allBars.length*timePoints) ) { // range = 500*timepoints\n";
		jScript << "				clearInterval(id);\n";
		jScript << "			} else {\n";
		jScript << "				counter = counter + speedup;\n";
		jScript << "				var newTime =  i_startTime + counter/(allBars.length*timePoints*(500/i_endTime));\n";
		jScript << "				document.getElementById('timeCount').innerHTML = (newTime).toFixed(2);\n";
		jScript << "				var cont = document.getElementById(\"tb\");\n";
		jScript << "				if (counter <= " << 500*d_barCount*i_bamFiles*turnoverTimes[1]/turnoverTimes[i_bamFiles-1] << ") {\n";
		jScript << "					pos[0][i] = pos[0][i] + speedup;\n";
		jScript << "					elem.style.height = baseLine[0][i] + pos[0][i]*posJumps[0][i] + 'px';\n";
		jScript << "					elem.style.top = 525 - baseLine[0][i] - pos[0][i]*posJumps[0][i] + 'px';\n"; 
		jScript << "				}\n";

		for (int k = 1; k < i_bamFiles-1; k++) {
			jScript << "				if ( (counter > " << 500*d_barCount*i_bamFiles*turnoverTimes[k]/turnoverTimes[i_bamFiles-1] << ") && (counter <= " << 500*d_barCount*i_bamFiles*turnoverTimes[k+1]/turnoverTimes[i_bamFiles-1] << ") ) {\n";
			jScript << "					pos[" << k << "][i] = pos[" << k << "][i] + speedup;\n";
			jScript << "					elem.style.height = baseLine[" << k << "][i] + pos[" << k << "][i]*posJumps[" << k << "][i] + 'px';\n";
			jScript << "					elem.style.top = 525 - baseLine[" << k << "][i] - pos[" << k << "][i]*posJumps[" << k << "][i] + 'px';\n"; 
			jScript << "				}\n";
		}

		jScript << "				if ( (counter > " << 500*d_barCount*i_bamFiles*turnoverTimes[i_bamFiles-1]/turnoverTimes[i_bamFiles-1] << ") && (counter <= " << 500*d_barCount*i_bamFiles*turnoverTimes[i_bamFiles]/turnoverTimes[i_bamFiles-1] << ") ) {\n";
		jScript << "					pos[" << i_bamFiles-1 << "][i] = pos[" << i_bamFiles-1 << "][i] + speedup;\n";
		jScript << "					elem.style.height = baseLine[" << i_bamFiles-1 << "][i] + pos[" << i_bamFiles-1 << "][i]*posJumps[" << i_bamFiles-1 << "][i] + 'px';\n";
		jScript << "					elem.style.top = 525 - baseLine[" << i_bamFiles-1 << "][i] - pos[" << i_bamFiles-1 << "][i]*posJumps[" << i_bamFiles-1 << "][i] + 'px';\n"; 
		jScript << "				}\n";
		jScript << "			}\n"; 
		jScript << "		}\n"; 
		jScript << "	}\n"; 
		jScript << "}\n"; 

		jScript << "function reset() {\n";
		jScript << "	for (i = 0; i < allBars.length; i++) {\n";
		jScript << "		var elem = document.getElementById(allBars[i]);\n";
		jScript << "		elem.style.top = 525-baseLine[0][i] + 'px';\n";
		jScript << "		elem.style.height = baseLine[0][i] + 'px';\n";
		jScript << "	}\n";
		jScript << "	document.getElementById(\"timeCount\").innerHTML = \"" << turnoverTimes[0] << "\"; // change to appropriate start time\n";
		jScript << "	time.style.left = 102 + 'px';\n";
		jScript << "}\n";

		jScript << "</script>\n";
		jScript << "</body>\n";
		jScript << "</html>\n";

		exonStartVector.clear();
		exonEndVector.clear();
		depthVector.clear();
		jScript.close();
	}
}
