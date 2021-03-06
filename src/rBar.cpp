#include <iostream>
#include <fstream>
#include <string>
using namespace std;

// creates a bar chart of absolute min and max values for all loci
void rBar(std::string s_rPltsName, int i_bamFiles, int turnoverTimes[], int maxIndexArray[], int minIndexArray[], int i_peakNumber) 
{

	// write R script
	ofstream rBar;
	rBar.open (s_rPltsName, std::ios::app);

	rBar << "time <- c(";
	for (int i = 0; i < i_bamFiles - 1; i++)
		rBar << "\"" << turnoverTimes[i] << "\", ";
	rBar <<  "\"" << turnoverTimes[i_bamFiles - 1] << "\")\n";

	// get a number to use as max Y axis value
	double d_yaxis = (double(minIndexArray[0]) / i_peakNumber * 100);

	rBar << "min <- c("; // min 
	for (int i = 0; i < i_bamFiles - 1; i++) {
		rBar << (double(minIndexArray[i]) / i_peakNumber * 100) << ", ";
		if ( (double(minIndexArray[i]) / i_peakNumber * 100) > d_yaxis) {
			d_yaxis = (double(minIndexArray[i]) / i_peakNumber * 100);
		}
	}
	rBar << (double(minIndexArray[i_bamFiles - 1]) / i_peakNumber * 100) << ")\n";
	if ( (double(minIndexArray[i_bamFiles - 1]) / i_peakNumber * 100) > d_yaxis) {
		d_yaxis = (double(minIndexArray[i_bamFiles - 1]) / i_peakNumber * 100);
	}

	rBar << "max <- c("; // max 
	for (int i = 0; i < i_bamFiles - 1; i++) {
		rBar << (double(maxIndexArray[i]) / i_peakNumber * 100) << ", ";
		if ( (double(maxIndexArray[i]) / i_peakNumber * 100) > d_yaxis) {
			d_yaxis = (double(maxIndexArray[i]) / i_peakNumber * 100);
		}
	}
	rBar << (double(maxIndexArray[i_bamFiles - 1]) / i_peakNumber * 100) << ")\n";
	if ( (double(maxIndexArray[i_bamFiles - 1]) / i_peakNumber * 100) > d_yaxis) {
		d_yaxis = (double(maxIndexArray[i_bamFiles - 1]) / i_peakNumber * 100);
	}


	rBar << "suppressMessages(dat1 <- melt(data.frame(min, max, time), variable.name=\"metric\"))\n";
	rBar << "# Do this so R doesn't alphabetically order x-axis\n";
	rBar << "dat1$time <- as.character(dat1$time)\n";
	rBar << "dat1$time <- factor(dat1$time, levels=unique(dat1$time))\n";
	rBar << "bar <- ggplot(dat1, aes(x=time, y=value, fill=metric)) +\n";
	rBar << "	geom_bar(stat=\"identity\", position=position_dodge(), colour=\"black\") +\n";
	rBar << "	scale_fill_manual(values=c(\"lemonchiffon3\", \"grey30\")) +\n";
	rBar << "	scale_y_continuous(limits = c(0, " << d_yaxis*1.15 << "), expand = c(0,0)) +\n";
	rBar << "	geom_text(data=dat1, aes(label=paste0(round(value,2),\"%\")), position=position_dodge(width=0.9), vjust=-0.25, size=2) +\n";
	rBar << "	labs(title = \"Occurence of normalized absolute minimum and maximum coverage\") +\n";
	rBar << "	labs(x = \"Time\",size = 10) +\n";
	rBar << "	labs(y = \"Percent Occurence\",size = 10) +\n";
	rBar << "	labs(fill = \"\") +\n";
	rBar << "	theme(\n";
	rBar << "		axis.text = element_text(size = 8, colour = \"black\"),\n";
	rBar << "		legend.key = element_rect(fill = \"white\"),\n";
	rBar << "		legend.background = element_rect(fill = \"white\"),\n";
	rBar << "		panel.grid.major = element_line(colour = \"white\"),\n";
	rBar << "		panel.grid.minor = element_blank(),\n";
	rBar << "		plot.title=element_text(size=14),\n";
	rBar << "		legend.key.size = unit(0.4, \"cm\"),\n";
	rBar << "		legend.text = element_text(size = 8),\n";
	rBar << "		panel.background = element_rect(fill = \"white\"),\n";
	rBar << "		axis.line.x = element_line(colour = \"black\"),\n";
	rBar << "		axis.line.y = element_line(colour = \"black\")\n";
	rBar << "	)\n";
	rBar.close();

}
