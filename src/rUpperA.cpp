#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h> 

// creates scatterplot boxplot of log2(rise/fall) of upper aymptotes in hills and valleys.
void rUpperA(std::string s_rPltsName, int i_peakNumber, int i_bamFiles, std::vector<std::vector<std::string> > &dataArray, int i_valley, int i_hill) 
{

/**	RECAP:
	dataArray[0][2*i_bamFiles + 3]  = "s_type_vec[i]"; 
	dataArray[0][2*i_bamFiles + 5]  = "s_drcFH_vec[i]"; 	// forward
	dataArray[0][2*i_bamFiles + 9]  = "s_drcRH_vec[i]"; 	// reverse

	dataArray[0][2*i_bamFiles + 6]  = "s_drcFU_vec[i]";	// forward
	dataArray[0][2*i_bamFiles + 10] = "s_drcRU_vec[i]";  	// reverse

	dataArray[i+1][2*i_bamFiles + 7]  = s_drcFL_vec[i];	// forward
	dataArray[i+1][2*i_bamFiles + 11] = s_drcRL_vec[i]; 	// reverse
**/

	int i_HV = i_hill + i_valley;

	// write R script
	std::ofstream rUpperA;
	rUpperA.open (s_rPltsName, std::ios::app);

	if (i_HV > 0) {

		// holds type: hills or valleys
		std::vector<std::string> s_upperTypeVec;
		std::vector<std::string> s_lowerTypeVec;
		// holds log2 ratio of rise/fall upper asymptotes
		std::vector<double> d_log2upperVec;
		std::vector<double> d_log2lowerVec;
	
		for (int i = 1; i < i_peakNumber+1; i++) {
			std::string s = dataArray[i][2*i_bamFiles + 3];
			if ( (s == "hill") || (s == "valley") )   { 
				double d_riseUpp = std::stod(dataArray[i][2*i_bamFiles + 6]);
				double d_fallUpp = std::stod(dataArray[i][2*i_bamFiles + 10]);
				if ( (d_riseUpp == 0) && (d_fallUpp == 0) ) {
					d_log2lowerVec.push_back(0);
					s_lowerTypeVec.push_back(s);
				} else if ( (isnan(log2(d_riseUpp/d_fallUpp))) || (isinf(log2(d_riseUpp/d_fallUpp))) ) {
				} else { 
					try {
						d_log2upperVec.push_back(log2(d_riseUpp/d_fallUpp)); 
						s_upperTypeVec.push_back(s);
					} catch(...){ }
				} 

			}
		}
		// LOWER ASYPTOTE
		for (int i = 1; i < i_peakNumber+1; i++) {
			std::string s = dataArray[i][2*i_bamFiles + 3];
			if ( (s == "hill") || (s == "valley") )   { 
				double d_riseUpp;
				double d_fallUpp;
				if (dataArray[i][2*i_bamFiles + 7] == "forced to 0") { d_riseUpp = 0; }
				else { d_riseUpp = std::stod(dataArray[i][2*i_bamFiles + 7]); }

				if (dataArray[i][2*i_bamFiles + 11] == "forced to 0") {  d_fallUpp = 0; }
				else { d_fallUpp = std::stod(dataArray[i][2*i_bamFiles + 11]); }

				if ( (d_riseUpp == 0) && (d_fallUpp == 0) ) { 
					d_log2lowerVec.push_back(0);
					s_lowerTypeVec.push_back(s);
				} else if ( (isnan(log2(d_riseUpp/d_fallUpp))) || (isinf(log2(d_riseUpp/d_fallUpp))) ) {
				} else {
					try {
						d_log2lowerVec.push_back(log2(d_riseUpp/d_fallUpp));
						s_lowerTypeVec.push_back(s);
					} catch(...){ }
				} 

			}
		}

		int i_iterator = s_upperTypeVec.size();
		if (i_iterator > 0) { 
			rUpperA << "log2upper <- c("; 
			for (int i = 0; i < i_iterator-1; i++) {
				rUpperA << d_log2upperVec[i] << ", ";	
			}		
			rUpperA << d_log2upperVec[i_iterator-1] << ")\n";	
	
			rUpperA << "logIDupper <- c("; 
			for (int i = 0; i < i_iterator-1; i++) {
				rUpperA << "\"" << s_upperTypeVec[i] << "\", ";	
			}		
			rUpperA << "\"" << s_upperTypeVec[i_iterator-1] << "\")\n";	
			rUpperA << "logUpDF <- data.frame(logIDupper,log2upper)\n";
			rUpperA << "UpperA <- ggplot(logUpDF , aes(logIDupper,log2upper, fill=logIDupper)) + \n";
			rUpperA << "	geom_boxplot(lwd=0.3,outlier.shape = NA) + geom_jitter(size = 0.1, alpha=0.15,width = 0.5) +\n";
			rUpperA << "	scale_fill_manual(values=c(\"white\", \"white\")) + \n";
			rUpperA << "	ggtitle(\"Upper Asymptote Symetricity\") +\n";
			rUpperA << "	labs(x = \"\") +\n";
			rUpperA << "	labs(y = \"log2(incline/decline)\",size = 8) +\n";
			rUpperA << "	theme(\n";
			rUpperA << "		legend.position=\"none\",\n";
			rUpperA << "		plot.title = element_text(size = 12),\n";
			rUpperA << "		axis.text = element_text(size = 8),\n";
			rUpperA << "		legend.key = element_rect(fill = \"white\"),\n";
			rUpperA << "		legend.background = element_rect(fill = \"white\"),\n";
			rUpperA << "		panel.grid.major = element_line(colour = \"white\"),\n";
			rUpperA << "		panel.grid.minor = element_blank(),\n";
			rUpperA << "		panel.background = element_rect(fill = \"white\"),\n";
			rUpperA << "		axis.line.x = element_line(colour = \"black\"),\n";
			rUpperA << "		axis.line.y = element_line(colour = \"black\")\n";
			rUpperA << "	)\n";
		} else {
			rUpperA << "#VARIABLE T\n";
			rUpperA << "d = data.frame(x1=c(1), x2=c(2), y1=c(1), y2=c(2), t=c('No hills or vallies'))\n";
			rUpperA << "#VARIABLE NAME\n";
			rUpperA << "UpperA <- ggplot() + \n";
			rUpperA << "scale_y_continuous(expand = c(0,0), limits = c(0, 2.5)) + \n";
			rUpperA << "scale_x_continuous(expand = c(0,0), limits = c(-2, 6)) + \n";
			rUpperA << "geom_text(data=d, aes(x=1.5, y=1.5, label=t), size=6) +\n";
			rUpperA << "	theme( \n";
			rUpperA << "		axis.text = element_text(size = 8, colour = \"white\"),\n";
			rUpperA << "		legend.key = element_rect(fill = \"white\"), \n";
			rUpperA << "		legend.background = element_rect(fill = \"white\"),\n";
			rUpperA << "		panel.grid.major = element_line(colour = \"white\"), \n";
			rUpperA << "		legend.position = \"\", \n";
			rUpperA << "		panel.grid.minor = element_blank(),\n";
			rUpperA << "		panel.background = element_rect(fill = \"white\"),\n";
			rUpperA << "		axis.line.x = element_line(colour = \"white\"),\n";
			rUpperA << "		axis.line.y = element_line(colour = \"white\"),\n";
			rUpperA << "		axis.title.x=element_blank(),\n";
			rUpperA << "		axis.text.x=element_blank(),\n";
			rUpperA << "		axis.ticks.x=element_blank(),\n";
			rUpperA << "		axis.title.y=element_blank(),\n";
			rUpperA << "		axis.text.y=element_blank(),\n";
			rUpperA << "		axis.ticks.y=element_blank()\n";
			rUpperA << "		) \n";
		}

		i_iterator = s_lowerTypeVec.size();
		if (i_iterator > 0) { 
			rUpperA << "log2lower <- c("; 
			for (int i = 0; i < i_iterator-1; i++) {
				rUpperA << d_log2lowerVec[i] << ", ";	
			}		
			rUpperA << d_log2lowerVec[i_iterator-1] << ")\n";

			rUpperA << "logIDlower <- c("; 
			for (int i = 0; i < i_iterator-1; i++) {
				rUpperA << "\"" << s_lowerTypeVec[i] << "\", ";	
			}		
			rUpperA << "\"" << s_lowerTypeVec[i_iterator-1] << "\")\n";	

			rUpperA << "logLowDF <- data.frame(logIDlower,log2lower)\n";
			rUpperA << "LowerA <- ggplot(logLowDF , aes(logIDlower,log2lower, fill=logIDlower)) + \n";
			rUpperA << "	geom_boxplot(lwd=0.3,outlier.shape = NA) + geom_jitter(size = 0.1, alpha=0.15,width = 0.5) +\n";
			rUpperA << "	scale_fill_manual(values=c(\"white\", \"white\")) + \n";
			rUpperA << "	ggtitle(\"Lower Asymptote Symetricity\") +\n";
			rUpperA << "	labs(x = \"\") +\n";
			rUpperA << "	labs(y = \"log2(incline/decline)\",size = 8) +\n";
			rUpperA << "	theme(\n";
			rUpperA << "		legend.position=\"none\",\n";
			rUpperA << "		plot.title = element_text(size = 12),\n";
			rUpperA << "		axis.text = element_text(size = 8),\n";
			rUpperA << "		legend.key = element_rect(fill = \"white\"),\n";
			rUpperA << "		legend.background = element_rect(fill = \"white\"),\n";
			rUpperA << "		panel.grid.major = element_line(colour = \"white\"),\n";
			rUpperA << "		panel.grid.minor = element_blank(),\n";
			rUpperA << "		panel.background = element_rect(fill = \"white\"),\n";
			rUpperA << "		axis.line.x = element_line(colour = \"black\"),\n";
			rUpperA << "		axis.line.y = element_line(colour = \"black\")\n";
			rUpperA << "	)\n";
		} else {
			rUpperA << "#VARIABLE T\n";
			rUpperA << "d = data.frame(x1=c(1), x2=c(2), y1=c(1), y2=c(2), t=c('No hills or vallies'))\n";
			rUpperA << "#VARIABLE NAME\n";
			rUpperA << "LowerA <- ggplot() + \n";
			rUpperA << "scale_y_continuous(expand = c(0,0), limits = c(0, 2.5)) + \n";
			rUpperA << "scale_x_continuous(expand = c(0,0), limits = c(-2, 6)) + \n";
			rUpperA << "geom_text(data=d, aes(x=1.5, y=1.5, label=t), size=6) +\n";
			rUpperA << "	theme( \n";
			rUpperA << "		axis.text = element_text(size = 8, colour = \"white\"),\n";
			rUpperA << "		legend.key = element_rect(fill = \"white\"), \n";
			rUpperA << "		legend.background = element_rect(fill = \"white\"),\n";
			rUpperA << "		panel.grid.major = element_line(colour = \"white\"), \n";
			rUpperA << "		legend.position = \"\", \n";
			rUpperA << "		panel.grid.minor = element_blank(),\n";
			rUpperA << "		panel.background = element_rect(fill = \"white\"),\n";
			rUpperA << "		axis.line.x = element_line(colour = \"white\"),\n";
			rUpperA << "		axis.line.y = element_line(colour = \"white\"),\n";
			rUpperA << "		axis.title.x=element_blank(),\n";
			rUpperA << "		axis.text.x=element_blank(),\n";
			rUpperA << "		axis.ticks.x=element_blank(),\n";
			rUpperA << "		axis.title.y=element_blank(),\n";
			rUpperA << "		axis.text.y=element_blank(),\n";
			rUpperA << "		axis.ticks.y=element_blank()\n";
			rUpperA << "		) \n";

		}

		s_upperTypeVec.clear();
		s_lowerTypeVec.clear();

	} else {
		rUpperA << "#VARIABLE T\n";
		rUpperA << "d = data.frame(x1=c(1), x2=c(2), y1=c(1), y2=c(2), t=c('No hills or vallies'))\n";
		rUpperA << "#VARIABLE NAME\n";
		rUpperA << "UpperA <- ggplot() + \n";
		rUpperA << "scale_y_continuous(expand = c(0,0), limits = c(0, 2.5)) + \n";
		rUpperA << "scale_x_continuous(expand = c(0,0), limits = c(-2, 6)) + \n";
		rUpperA << "geom_text(data=d, aes(x=1.5, y=1.5, label=t), size=6) +\n";
		rUpperA << "	theme( \n";
		rUpperA << "		axis.text = element_text(size = 8, colour = \"white\"),\n";
		rUpperA << "		legend.key = element_rect(fill = \"white\"), \n";
		rUpperA << "		legend.background = element_rect(fill = \"white\"),\n";
		rUpperA << "		panel.grid.major = element_line(colour = \"white\"), \n";
		rUpperA << "		legend.position = \"\", \n";
		rUpperA << "		panel.grid.minor = element_blank(),\n";
		rUpperA << "		panel.background = element_rect(fill = \"white\"),\n";
		rUpperA << "		axis.line.x = element_line(colour = \"white\"),\n";
		rUpperA << "		axis.line.y = element_line(colour = \"white\"),\n";
		rUpperA << "		axis.title.x=element_blank(),\n";
		rUpperA << "		axis.text.x=element_blank(),\n";
		rUpperA << "		axis.ticks.x=element_blank(),\n";
		rUpperA << "		axis.title.y=element_blank(),\n";
		rUpperA << "		axis.text.y=element_blank(),\n";
		rUpperA << "		axis.ticks.y=element_blank()\n";
		rUpperA << "		) \n";

		rUpperA << "#VARIABLE T\n";
		rUpperA << "d = data.frame(x1=c(1), x2=c(2), y1=c(1), y2=c(2), t=c('No hills or vallies'))\n";
		rUpperA << "#VARIABLE NAME\n";
		rUpperA << "LowerA <- ggplot() + \n";
		rUpperA << "scale_y_continuous(expand = c(0,0), limits = c(0, 2.5)) + \n";
		rUpperA << "scale_x_continuous(expand = c(0,0), limits = c(-2, 6)) + \n";
		rUpperA << "geom_text(data=d, aes(x=1.5, y=1.5, label=t), size=6) +\n";
		rUpperA << "	theme( \n";
		rUpperA << "		axis.text = element_text(size = 8, colour = \"white\"),\n";
		rUpperA << "		legend.key = element_rect(fill = \"white\"), \n";
		rUpperA << "		legend.background = element_rect(fill = \"white\"),\n";
		rUpperA << "		panel.grid.major = element_line(colour = \"white\"), \n";
		rUpperA << "		legend.position = \"\", \n";
		rUpperA << "		panel.grid.minor = element_blank(),\n";
		rUpperA << "		panel.background = element_rect(fill = \"white\"),\n";
		rUpperA << "		axis.line.x = element_line(colour = \"white\"),\n";
		rUpperA << "		axis.line.y = element_line(colour = \"white\"),\n";
		rUpperA << "		axis.title.x=element_blank(),\n";
		rUpperA << "		axis.text.x=element_blank(),\n";
		rUpperA << "		axis.ticks.x=element_blank(),\n";
		rUpperA << "		axis.title.y=element_blank(),\n";
		rUpperA << "		axis.text.y=element_blank(),\n";
		rUpperA << "		axis.ticks.y=element_blank()\n";
		rUpperA << "		) \n";
	}

	

	rUpperA.close();

}
