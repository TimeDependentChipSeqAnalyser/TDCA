#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

// creates scatter plots of inflection points and hills coef and binned density plots of inf and hills coef 
void rStdev(std::string s_rPltsName, int i_bamFiles, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, int turnoverTimes[], int i_undefRise, int i_undefFall, int i_valley, int i_hill, int i_fall, int i_rise) 
{

	// write R script
	ofstream rStdev;
	rStdev.open (s_rPltsName, std::ios::app);

	// write csv files // XXX ADDED XXX
	ofstream FcsvInf;
	ofstream FcsvStdev;
	ofstream RcsvInf;
	ofstream RcsvStdev;
	ofstream stdevBoxPlot;

	std::string s_FcsvInfName;
	std::string s_FcsvStdevName;
	std::string s_RcsvInfName;
	std::string s_RcsvStdevName;
	std::string s_stdevBoxPlotName;

	s_FcsvInfName = "./" + s_rPltsName + "_data/fRepInf.csv"; 
	s_RcsvInfName = "./" + s_rPltsName + "_data/rRepInf.csv"; 
	s_FcsvStdevName = "./" + s_rPltsName + "_data/fRepStdev.csv"; 
	s_RcsvStdevName = "./" + s_rPltsName + "_data/rRepStdev.csv";
	s_stdevBoxPlotName = "./" + s_rPltsName + "_data/stdevBox.csv";

	FcsvInf.open(s_FcsvInfName);
	RcsvInf.open(s_RcsvInfName);
	FcsvStdev.open(s_FcsvStdevName);
	RcsvStdev.open(s_RcsvStdevName);
	stdevBoxPlot.open(s_stdevBoxPlotName);

	int i_increase = i_rise + i_undefRise + i_hill + i_valley;
	int i_decrease = i_fall + i_undefFall + i_hill + i_valley;
	int i_iterator = 0;

	double d_maxFstdev = 0;
	double d_maxRstdev = 0;

	if (i_increase > 0) {
		// increase
		for (int i = 0; i < i_peakNumber; i++) {
			std::string s_line = dataArray[i+1][2*i_bamFiles + 5];
			if ( (s_line != "-") &&  (s_line != "NaN") ) {
				if (i_iterator == i_increase-1) {  // last data value
					FcsvInf << s_line;
					i_iterator++;
				} 
				if (i_iterator < i_increase-1)  {  // not last data value
					FcsvInf << s_line << "\n";
					i_iterator++;
				}  
			}
		}

		i_iterator = 0;
		for (int i = 0; i < i_peakNumber; i++) {
			std::string s_line = dataArray[i+1][2*i_bamFiles + 5]; 
			if ( (s_line != "-") &&  (s_line != "NaN") ) {
				double d_sum = 0;
				if (i_iterator == i_increase-1) { // last data value
					for (int k = 0; k < i_bamFiles; k++) { 	// for each time point
						d_sum += stod(dataArray[i+1][3 + i_bamFiles + k]);
					}
					FcsvStdev << d_sum/i_bamFiles; 
					i_iterator++;
					if ( (d_sum/i_bamFiles) > d_maxFstdev) { d_maxFstdev = (d_sum/i_bamFiles); }
				}
				if (i_iterator < i_increase-1) { // not last data value
					for (int k = 0; k < i_bamFiles; k++) { 	// for each time point
						d_sum += stod(dataArray[i+1][3 + i_bamFiles + k]);
					}
					FcsvStdev << d_sum/i_bamFiles << "\n"; 
					i_iterator++;
					if ( (d_sum/i_bamFiles) > d_maxFstdev) { d_maxFstdev = (d_sum/i_bamFiles); }
				}			
			}
		}
	}

	i_iterator = 0;
	if (i_decrease > 0) {
		// increase
		for (int i = 0; i < i_peakNumber; i++) {
			std::string s_line = dataArray[i+1][2*i_bamFiles + 12];
			if ( (s_line != "-") &&  (s_line != "NaN") ) {
				
				if (i_iterator == i_decrease-1) { // last data value
					RcsvInf << s_line; 
					i_iterator++;
				} 
				if (i_iterator < i_decrease-1)  { // not last data value
					RcsvInf << s_line << "\n"; 
					i_iterator++;
				}  
			}
		}

		i_iterator = 0;
		for (int i = 0; i < i_peakNumber; i++) {
			std::string s_line = dataArray[i+1][2*i_bamFiles + 12]; 
			if ( (s_line != "-") &&  (s_line != "NaN") ) {
				double d_sum = 0;
				if (i_iterator == i_decrease-1) { // last data value
					for (int k = 0; k < i_bamFiles; k++) { 	// for each time point
						d_sum += stod(dataArray[i+1][3 + i_bamFiles + k]);
					}
					RcsvStdev << d_sum/i_bamFiles; 
					i_iterator++;
					if ( (d_sum/i_bamFiles) > d_maxRstdev) { d_maxRstdev = (d_sum/i_bamFiles); }
				}
				if (i_iterator < i_decrease-1) { // not last data value
					for (int k = 0; k < i_bamFiles; k++) { 	// for each time point
						d_sum += stod(dataArray[i+1][3 + i_bamFiles + k]);
					}
					RcsvStdev << d_sum/i_bamFiles << "\n"; 
					i_iterator++;
					if ( (d_sum/i_bamFiles) > d_maxRstdev) { d_maxRstdev = (d_sum/i_bamFiles); }
				}			
			}
		}
	}

	if (i_increase > 0) {
		rStdev << "fInfRep <- read.csv(file=\"./" + s_rPltsName + "_data/fRepInf.csv\", header=FALSE)\n";
		rStdev << "fStdevRep <- read.csv(file=\"./" + s_rPltsName + "_data/fRepStdev.csv\", header=FALSE)\n";
		rStdev << "fStdev<- data.frame(fStdevRep ,fInfRep)\n"; 
		rStdev << "s1 <- ggplot(fStdev, aes(fStdev[,2],fStdev[,1])) +\n"; 
		rStdev << "	geom_point(size = 0.3) +\n"; 
		rStdev << "	ggtitle(\"Coverage Standard Deviation and TTI of Signal Increase Loci\") +\n"; 
		rStdev << "	labs(x = \"TTI (time)\",size = 8) +\n"; 
		rStdev << "	labs(y = \"Coverage Standard Deviation (reads)\",size = 8) +\n"; 
		rStdev << "	scale_y_continuous(expand = c(0,0), limits = c(0, " << d_maxFstdev*0.9 << ")) +\n";  
		rStdev << "	scale_x_continuous(expand = c(0,0), limits = c(0, " << turnoverTimes[i_bamFiles-1] << ")) +\n";  
		rStdev << "	theme(\n"; 
		rStdev << "		axis.title.x=element_text(size=8),\n"; 
		rStdev << "		axis.title.y=element_text(size=8),\n"; 
		rStdev << "		plot.title=element_text(size=10),\n"; 
		rStdev << "		axis.text = element_text(size = 8, colour = \"black\"),\n"; 
		rStdev << "		legend.key = element_rect(fill = \"white\"),\n"; 
		rStdev << "		legend.background = element_rect(fill = \"white\"),\n"; 
		rStdev << "		panel.grid.major = element_line(colour = \"white\"),\n"; 
		rStdev << "		panel.grid.minor = element_blank(),\n"; 
		rStdev << "		panel.background = element_rect(fill = \"white\"),\n"; 
		rStdev << "		axis.line.x = element_line(colour = \"black\"),\n"; 
		rStdev << "		axis.line.y = element_line(colour = \"black\")\n"; 
		rStdev << "	)\n"; 
	} else {
		rStdev << "d = data.frame(x1=c(1), x2=c(2), y1=c(1), y2=c(2), t=c('No rises'))\n";
		rStdev << "s1 <- ggplot() + \n";
		rStdev << "scale_y_continuous(expand = c(0,0), limits = c(0, 2.5)) + \n";
		rStdev << "scale_x_continuous(expand = c(0,0), limits = c(-2, 6)) + \n";
		rStdev << "geom_text(data=d, aes(x=1.5, y=1.5, label=t), size=6) +\n";
		rStdev << "	theme( \n";
		rStdev << "		axis.text = element_text(size = 8, colour = \"white\"),\n";
		rStdev << "		legend.key = element_rect(fill = \"white\"), \n";
		rStdev << "		legend.background = element_rect(fill = \"white\"),\n";
		rStdev << "		panel.grid.major = element_line(colour = \"white\"), \n";
		rStdev << "		legend.position = \"\", \n";
		rStdev << "		panel.grid.minor = element_blank(),\n";
		rStdev << "		panel.background = element_rect(fill = \"white\"),\n";
		rStdev << "		axis.line.x = element_line(colour = \"white\"),\n";
		rStdev << "		axis.line.y = element_line(colour = \"white\"),\n";
		rStdev << "		axis.title.x=element_blank(),\n";
		rStdev << "		axis.text.x=element_blank(),\n";
		rStdev << "		axis.ticks.x=element_blank(),\n";
		rStdev << "		axis.title.y=element_blank(),\n";
		rStdev << "		axis.text.y=element_blank(),\n";
		rStdev << "		axis.ticks.y=element_blank()\n";
		rStdev << "		) \n";
	}

	if (i_decrease > 0) {
		rStdev << "rInfRep <- read.csv(file=\"./" + s_rPltsName + "_data/rRepInf.csv\", header=FALSE)\n";
		rStdev << "rStdevRep <- read.csv(file=\"./" + s_rPltsName + "_data/rRepStdev.csv\", header=FALSE)\n";
		rStdev << "rStdev<- data.frame(rStdevRep,rInfRep)\n"; 
		rStdev << "s2 <- ggplot(rStdev, aes(rStdev[,2],rStdev[,1])) +\n"; 
		rStdev << "	geom_point(size = 0.3) +\n"; 
		rStdev << "	ggtitle(\"Coverage Standard Deviation and TTI of Signal Decrease Loci\") +\n"; 
		rStdev << "	labs(x = \"TTI (time)\",size = 8) +\n"; 
		rStdev << "	labs(y = \"Coverage Standard Deviation (reads)\",size = 8) +\n"; 
		rStdev << "	scale_y_continuous(expand = c(0,0), limits = c(0, " << d_maxRstdev*0.9 << ")) +\n"; 
		rStdev << "	scale_x_continuous(expand = c(0,0), limits = c(0, " << turnoverTimes[i_bamFiles-1] << ")) +\n"; 
		rStdev << "	theme(\n"; 
		rStdev << "		axis.title.x=element_text(size=8),\n"; 
		rStdev << "		axis.title.y=element_text(size=8),\n"; 
		rStdev << "		plot.title=element_text(size=10),\n"; 
		rStdev << "		axis.text = element_text(size = 8, colour = \"black\"),\n"; 
		rStdev << "		legend.key = element_rect(fill = \"white\"),\n"; 
		rStdev << "		legend.background = element_rect(fill = \"white\"),\n"; 
		rStdev << "		panel.grid.major = element_line(colour = \"white\"),\n"; 
		rStdev << "		panel.grid.minor = element_blank(),\n"; 
		rStdev << "		panel.background = element_rect(fill = \"white\"),\n"; 
		rStdev << "		axis.line.x = element_line(colour = \"black\"),\n"; 
		rStdev << "		axis.line.y = element_line(colour = \"black\")\n"; 
		rStdev << "	)\n"; 
	} else {
		rStdev << "d = data.frame(x1=c(1), x2=c(2), y1=c(1), y2=c(2), t=c('No falls'))\n";
		rStdev << "s2 <- ggplot() + \n";
		rStdev << "scale_y_continuous(expand = c(0,0), limits = c(0, 2.5)) + \n";
		rStdev << "scale_x_continuous(expand = c(0,0), limits = c(-2, 6)) + \n";
		rStdev << "geom_text(data=d, aes(x=1.5, y=1.5, label=t), size=6) +\n";
		rStdev << "	theme( \n";
		rStdev << "		axis.text = element_text(size = 8, colour = \"white\"),\n";
		rStdev << "		legend.key = element_rect(fill = \"white\"), \n";
		rStdev << "		legend.background = element_rect(fill = \"white\"),\n";
		rStdev << "		panel.grid.major = element_line(colour = \"white\"), \n";
		rStdev << "		legend.position = \"\", \n";
		rStdev << "		panel.grid.minor = element_blank(),\n";
		rStdev << "		panel.background = element_rect(fill = \"white\"),\n";
		rStdev << "		axis.line.x = element_line(colour = \"white\"),\n";
		rStdev << "		axis.line.y = element_line(colour = \"white\"),\n";
		rStdev << "		axis.title.x=element_blank(),\n";
		rStdev << "		axis.text.x=element_blank(),\n";
		rStdev << "		axis.ticks.x=element_blank(),\n";
		rStdev << "		axis.title.y=element_blank(),\n";
		rStdev << "		axis.text.y=element_blank(),\n";
		rStdev << "		axis.ticks.y=element_blank()\n";
		rStdev << "		) \n";
	}


	// write data to stdevBoxPlot
	for (int i = 1; i < i_peakNumber; i++) { 		// for each peak
		for (int k = 0; k < i_bamFiles; k++) { 	// for each time point
			stdevBoxPlot << dataArray[i][3 + i_bamFiles + k] << "\n";
		}
	}
	for (int k = 0; k < i_bamFiles-1; k++)  		// last peak
		stdevBoxPlot << dataArray[i_peakNumber][3 + i_bamFiles + k] << "\n";
	stdevBoxPlot << dataArray[i_peakNumber][3 + i_bamFiles + i_bamFiles-1];


	rStdev << "stdevvalue <- read.csv(file=\"./" + s_rPltsName + "_data/stdevBox.csv\", header=FALSE)\n";

	rStdev << "stdevtime <- c (";
	for (int k = 0; k < i_bamFiles-1; k++)  	
		rStdev << "\"" << turnoverTimes[k] << "\", ";
	rStdev << "\"" << turnoverTimes[i_bamFiles-1] << "\")\n";

	rStdev << "stdev <- data.frame(stdevvalue, stdevtime)\n";
	rStdev << "# Do this so R doesn't alphabetically order x-axis\n";
	rStdev << "stdev$stdevtime <- as.character(stdev$stdevtime)\n";
	rStdev << "stdev$stdevtime <- factor(stdev$stdevtime, levels=unique(stdev$stdevtime))\n";
	rStdev << "stdevBox <- ggplot(stdev, aes(x=stdev[,2], y=stdev[,1], fill=stdev[,2])) +\n";
	rStdev << "	geom_boxplot(outlier.shape = NA) +\n";
	rStdev << "	scale_fill_brewer(palette=\"GnBu\") +\n";
	rStdev << "	ggtitle(\"Coverage Standard Deviation of Time Points\") +\n";
	rStdev << "	labs(y = \"Coverage Standard Deviation (reads)\",size = 8) +\n";
	rStdev << "	labs(x = \"Time Points\",size = 8) +\n";
	rStdev << "	coord_cartesian(ylim = range(boxplot(stdev[,1], plot=FALSE)$stats)*c(.9, 2)) +\n";
	rStdev << "	theme(\n";
	rStdev << "		axis.title.x=element_text(size=8),\n"; 
	rStdev << "		axis.title.y=element_text(size=8),\n"; 
	rStdev << "		plot.title=element_text(size=10),\n"; 
	rStdev << "		axis.text = element_text(size = 8, colour = \"black\"),\n";
	rStdev << "		legend.position=(\"none\"),\n";
	rStdev << "		panel.grid.major = element_line(colour = \"white\"),\n";
	rStdev << "		panel.grid.minor = element_blank(),\n";
	rStdev << "		panel.background = element_rect(fill = \"white\"),\n";
	rStdev << "		axis.line.x = element_line(colour = \"black\"),\n";
	rStdev << "		axis.line.y = element_line(colour = \"black\")\n";
	rStdev << "	)\n";

	rStdev.close();
	FcsvInf.close();
	FcsvStdev.close();
	RcsvInf.close();
	RcsvStdev.close();
	stdevBoxPlot.close();
	//csvStdevTime.close();
	//csvStdevALL.close();
}
