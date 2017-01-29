#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <math.h>
#include <array>

std::string exec(const char* cmd);

void rHeatmap(std::string s_rPltsName, int i_peakNumber, int i_bamFiles, std::vector<std::vector<std::string> > &dataArray, int turnoverTimes[]) {

	// output
	std::ofstream rHeatmap;
	rHeatmap.open (s_rPltsName, std::ios::app);

	std::string s_csvX1Name = "./" + s_rPltsName + "_data/heatmapX1.csv";
	std::string s_csvX2Name = "./" + s_rPltsName + "_data/heatmapX2.csv";
	std::string s_csvY1Name = "./" + s_rPltsName + "_data/heatmapY1.csv";
	std::string s_csvY2Name = "./" + s_rPltsName + "_data/heatmapY2.csv";
	std::string s_csvDepthName = "./" + s_rPltsName + "_data/heatmapDepth.csv";

	std::ofstream csvX1;
	std::ofstream csvX2;
	std::ofstream csvY1;
	std::ofstream csvY2;
	std::ofstream csvDepth;

	csvX1.open(s_csvX1Name);
	csvX2.open(s_csvX2Name);
	csvY1.open(s_csvY1Name);
	csvY2.open(s_csvY2Name);
	csvDepth.open(s_csvDepthName);



	// print normalized depth values to CSV file
	double d_max;
	double d_min;
	double d_depthDiff;
	for (int i = 1; i < i_peakNumber; i++) {
		d_max = stod(dataArray[i][3]);
		d_min = stod(dataArray[i][3]);
		for (int j = 0; j < i_bamFiles; j++) { // establish min and max at each loci
			if (stod(dataArray[i][j+3]) > d_max) { d_max = stod(dataArray[i][j+3]); }
			if (stod(dataArray[i][j+3]) < d_min) { d_min = stod(dataArray[i][j+3]); }
		}
		d_depthDiff = d_max - d_min;
		for (int j = 0; j < i_bamFiles; j++) {
			if (stod(dataArray[i][j+3]) - d_min == 0) { csvDepth << 0 << '\n'; }
			else { csvDepth << (stod(dataArray[i][j+3]) - d_min ) / d_depthDiff << '\n'; }
		}		
	}
	// last peak
	d_max = stod(dataArray[i_peakNumber][3]);
	d_min = stod(dataArray[i_peakNumber][3]);
	for (int j = 0; j < i_bamFiles; j++) {
		if (stod(dataArray[i_peakNumber][j+3]) > d_max) { d_max = stod(dataArray[i_peakNumber][j+3]); }
		if (stod(dataArray[i_peakNumber][j+3]) < d_min) { d_min = stod(dataArray[i_peakNumber][j+3]); }
	}
	d_depthDiff = d_max - d_min;
	for (int j = 0; j < (i_bamFiles-1); j++) {
		if (stod(dataArray[i_peakNumber][j+3]) - d_min == 0) { csvDepth << 0 << '\n'; }
		else { csvDepth << (stod(dataArray[i_peakNumber][j+3]) - d_min ) / d_depthDiff << '\n'; }
	}
	if (stod(dataArray[i_peakNumber][i_bamFiles+2]) - d_min == 0) { csvDepth << 0; }
	else { csvDepth << (stod(dataArray[i_peakNumber][i_bamFiles+2]) - d_min ) / d_depthDiff; }




	// print X1 and X2 to CSV
	for (int i = 0; i < i_peakNumber-1; i++) {	
		for (int j = 1; j < i_bamFiles+1; j++) {
			csvX1 << j << '\n'; 
			csvX2 << j+1 << '\n'; 
		}
	}
	// last peak
	for (int j = 1; j < i_bamFiles; j++) {
		csvX1 << j << '\n'; 
		csvX2 << j+1 << '\n'; 
	}
	csvX1 << i_bamFiles;
	csvX2 << i_bamFiles+1;

	// print Y1 and Y2 to CSV
	for (int i = i_peakNumber-1; i > 0; i--) {	
		for (int j = 1; j < i_bamFiles+1; j++) {
			csvY1 << i+0.05 << '\n'; 
			csvY2 << i+0.95 << '\n'; 
		}
	}
	// last peak
	for (int j = 1; j < i_bamFiles; j++) {
		csvY1 << 0.05 << '\n'; 
		csvY2 << 0.95 << '\n'; 
	}
	csvY1 << 1.05; 
	csvY2 << 1.95; 



	rHeatmap << "heatmapX1 <- read.csv(file=\"./" << s_rPltsName << "_data/heatmapX1.csv\", header=FALSE)\n";
	rHeatmap << "heatmapX2 <- read.csv(file=\"./" << s_rPltsName << "_data/heatmapX2.csv\", header=FALSE)\n";
	rHeatmap << "heatmapY1 <- read.csv(file=\"./" << s_rPltsName << "_data/heatmapY1.csv\", header=FALSE)\n";
	rHeatmap << "heatmapY2 <- read.csv(file=\"./" << s_rPltsName << "_data/heatmapY2.csv\", header=FALSE)\n";
	rHeatmap << "heatmapFill <- read.csv(file=\"./" << s_rPltsName << "_data/heatmapDepth.csv\", header=FALSE)\n";
	rHeatmap << "heatmapFrame <- data.frame(heatmapX1, heatmapX2, heatmapY1, heatmapY2, heatmapFill)\n";

	rHeatmap << "heatmap = ggplot() + \n";
	rHeatmap << "scale_x_continuous(limits = c(0, " << i_bamFiles+2 << "), name=\"x\") + \n";
	rHeatmap << "scale_y_continuous(limits = c(0, " << int(i_peakNumber*1.1) << ")) +\n";
	rHeatmap << "geom_rect(data=heatmapFrame, mapping=aes(xmin=heatmapFrame[,1], xmax=heatmapFrame[,2], ymin=heatmapFrame[,3], ymax=heatmapFrame[,4], fill=heatmapFrame[,5])) +\n";
	rHeatmap << "scale_fill_gradient(low=\"white\", high=\"darkblue\") +\n";
	// box
	rHeatmap << "geom_rect(mapping=aes(xmin=1, xmax=" << i_bamFiles+1 << ", ymin=0, ymax=" << i_peakNumber << "), color=\"black\", alpha=0) +\n";
	// time labels
	for (int j = 0; j < i_bamFiles; j++) {
		rHeatmap << "geom_text(aes(x=" << j+1 << ".5, y=" << i_peakNumber*1.02 << ", label=paste(\"" << turnoverTimes[j] << "\"))) +\n";
	}
	rHeatmap << "ggtitle(\"Normalized depth heatmap at each loci\") +\n";
	// lines
	//rHeatmap << "geom_rect(mapping=aes(xmin=2, xmax=2, ymin=0, ymax=" << i_wc << "), color=\"black\", alpha=0) +\n";
	//rHeatmap << "geom_rect(mapping=aes(xmin=3, xmax=3, ymin=0, ymax=" << i_wc << "), color=\"black\", alpha=0) +\n";
	//rHeatmap << "geom_rect(mapping=aes(xmin=4, xmax=4, ymin=0, ymax=" << i_wc << "), color=\"black\", alpha=0) +\n";
	// scale
	rHeatmap << "geom_rect(mapping=aes(xmin=0.5, xmax=0.5, ymin=0, ymax=" << int(i_peakNumber*0.1) << "), color=\"black\", alpha=0) +\n";
	rHeatmap << "geom_text(aes(x=0.3, y=" << int(i_peakNumber*0.1)/2 << ", label=paste(\"" << int(i_peakNumber*0.1) << " loci\"), angle = 90)) +\n";

	rHeatmap << "geom_rect(mapping=aes(xmin=0.25, xmax=0.75, ymin=" << int(i_peakNumber*0.3) << ", ymax=" << int(i_peakNumber*0.4) << "), color=\"black\", fill=\"white\") +\n";
	rHeatmap << "geom_text(aes(x=0.05, y=" << int(i_peakNumber*0.35) << ", label=paste(\"low\"), angle = 90)) +\n";

	rHeatmap << "geom_rect(mapping=aes(xmin=0.25, xmax=0.75, ymin=" << int(i_peakNumber*0.4) << ", ymax=" << int(i_peakNumber*0.5) << "), color=\"black\", fill = \"darkblue\") +\n";
	rHeatmap << "geom_text(aes(x=0.05, y=" << int(i_peakNumber*0.45) << ", label=paste(\"high\"), angle = 90)) +\n";


	rHeatmap << "theme(legend.position=\"none\",\n";
	rHeatmap << "legend.background = element_rect(fill = \"white\"),\n";
	rHeatmap << "panel.grid.major = element_line(colour = \"white\"),\n";
	rHeatmap << "panel.grid.minor = element_blank(),\n";
	rHeatmap << "panel.background = element_rect(fill = \"white\"),\n";
	rHeatmap << "axis.line.x = element_line(colour = \"white\"),\n";
	rHeatmap << "axis.title.x=element_blank(),\n";
	rHeatmap << "axis.text.x=element_blank(),\n";
	rHeatmap << "axis.ticks.x=element_blank(),\n";
	rHeatmap << "axis.title.y=element_blank(),\n";
	rHeatmap << "axis.text.y=element_blank(),\n";
	rHeatmap << "axis.ticks.y=element_blank(),\n";
	rHeatmap << "axis.line.y = element_line(colour = \"white\")) \n\n";




	csvX1.close();
	csvX2.close();
	csvY1.close();
	csvY2.close();
	csvDepth.close();
	rHeatmap.close();
}










