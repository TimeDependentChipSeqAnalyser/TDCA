#include <iostream>
#include <fstream>
#include <string>
#include <vector> 

// creates normalize averaged profiles
void rProfiles(std::string s_rPltsName, int i_peakNumber, int i_bamFiles, std::vector<std::vector<std::string> > &dataArray, bool b_model, double d_satthresh, int turnoverTimes[], int i_valley, int i_hill, int i_fall, int i_rise, int i_undefRise, int i_undefFall) 
{

/**	RECAP:
	dataArray[0][2*i_bamFiles + 3]  = "s_type_vec[i]"; 
	dataArray[0][2*i_bamFiles + 5]  = "s_drcFH_vec[i]";
	dataArray[0][2*i_bamFiles + 6]  = "s_drcFU_vec[i]"; 
	dataArray[0][2*i_bamFiles + 9]  = "s_drcRH_vec[i]"; 
	dataArray[0][2*i_bamFiles + 10] = "s_drcRU_vec[i]"; 
**/

	// write R script
	std::ofstream rProfiles;
	rProfiles.open (s_rPltsName, std::ios::app);

	// [rows][columns]
	double d_riseArr[i_bamFiles];
	double d_hillArr[i_bamFiles];
	double d_valleyArr[i_bamFiles];
	double d_fallArr[i_bamFiles];
	double d_undefRiseArr[i_bamFiles];
	double d_undefFallArr[i_bamFiles];

	for (int i = 0; i < i_bamFiles; i++) {
		d_riseArr[i]      = 0;
		d_hillArr[i]      = 0;
		d_valleyArr[i]    = 0;
		d_fallArr[i]      = 0;
		d_undefRiseArr[i] = 0;
		d_undefFallArr[i] = 0;
	}

	int i_riseCount = 0;
	int i_hillCount = 0;
	int i_valleyCount = 0;
	int i_fallCount = 0;
	int i_undefRiseCount = 0;
	int i_undefFallCount = 0;

	for (int i = 1; i < i_peakNumber+1; i++) {
		double d_max = stod(dataArray[i][3]);
		double d_min = stod(dataArray[i][3]);
		for (int j = 1; j < i_bamFiles; j++) { // establish min/max
			if (stod(dataArray[i][j+3]) > d_max) { d_max = stod(dataArray[i][j+3]); } 
			if (stod(dataArray[i][j+3]) < d_min) { d_min = stod(dataArray[i][j+3]); } 
		}
		std::string s_t = dataArray[i][2*i_bamFiles + 3];
		if (s_t == "rise") {
			for (int j = 0; j < i_bamFiles; j++) {
				d_riseArr[j] += (stod(dataArray[i][j+3])-d_min)/d_max*100;
				i_riseCount++;
			}
		}
		if (s_t == "hill") {
			for (int j = 0; j < i_bamFiles; j++) {
				d_hillArr[j] += (stod(dataArray[i][j+3])-d_min)/d_max*100;
				i_hillCount++;
			}
		}
		if (s_t == "valley") {
			for (int j = 0; j < i_bamFiles; j++) {
				d_valleyArr[j] += (stod(dataArray[i][j+3])-d_min)/d_max*100;
				i_valleyCount++;
			}
		}
		if (s_t == "fall") {
			for (int j = 0; j < i_bamFiles; j++) {
				d_fallArr[j] += (stod(dataArray[i][j+3])-d_min)/d_max*100;
				i_fallCount++;
			}
		}
		std::string s_FHC = dataArray[i][2*i_bamFiles + 5];
		if ( (s_t == "undefined") && (s_FHC != "NaN") && (s_FHC != "-") ) { // undefined rise
			for (int j = 0; j < i_bamFiles; j++) {
				d_undefRiseArr[j] += (stod(dataArray[i][j+3])-d_min)/d_max*100;
				i_undefRiseCount++;
			}
		}
		std::string s_RHC = dataArray[i][2*i_bamFiles + 9];
		if ( (s_t == "undefined") && (s_RHC != "NaN") && (s_RHC != "-") ) { // undefined fall
			for (int j = 0; j < i_bamFiles; j++) {
				d_undefFallArr[j] += (stod(dataArray[i][j+3])-d_min)/d_max*100;
				i_undefFallCount++;
			}
		}	
	}

	double d_thresh = 1-d_satthresh;

	// hills
	double d_max = d_hillArr[0];
	double d_min = d_hillArr[0];
	int i_maxIndex = 0;
	for (int j = 1; j < i_bamFiles; j++) { // establish min/max
		if (d_hillArr[j] > d_max) { d_max = d_hillArr[j]; i_maxIndex = j; } 
		if (d_hillArr[j] < d_min) { d_min = d_hillArr[j]; } 
	}

	int i_maxTrail = 0;
	for (int j = i_maxIndex - 1; j > -1; j--) {
		if ( d_hillArr[j] >= (d_max-((d_max-d_min)*d_thresh)) ) { i_maxTrail++; }	
	}

	int i_maxLead = 0;
	for (int j = i_maxIndex + 1; j < i_bamFiles; j++) {
		if ( d_hillArr[j] >= (d_max-((d_max-d_min)*d_thresh)) ) { i_maxLead++; }
	}

	// valleys
	d_max = d_valleyArr[0];
	d_min = d_valleyArr[0];
	int i_minIndex = 0;
	for (int j = 1; j < i_bamFiles; j++) { // establish min/max
		if (d_valleyArr[j] > d_max) { d_max = d_valleyArr[j]; } 
		if (d_valleyArr[j] < d_min) { d_min = d_valleyArr[j]; i_minIndex = j; } 
	}

	int i_minTrail = 0;
	for (int j = i_minIndex - 1; j > -1; j--) {
		if ( d_valleyArr[j] <= ( d_min+((d_max-d_min)*d_thresh)) ) { i_minTrail++; }		
	}

	int i_minLead = 0;
	for (int j = i_minIndex + 1; j < i_bamFiles; j++) {
		if ( d_valleyArr[j] <= ( d_min+((d_max-d_min)*d_thresh)) ) { i_minLead++; }		
	}

	// Print Values
	if (b_model) {
		if (i_hill > 0) {
			// HILLS
			rProfiles << "hillRiseTime <- c(";
			for (int i = 0; i < i_maxIndex+i_maxLead; i++) {
				rProfiles << turnoverTimes[i] << ", ";
			}
			rProfiles << turnoverTimes[i_maxIndex+i_maxLead] << ")\n";

			rProfiles << "hillRiseDepth <- c(";
			for (int i = 0; i < i_maxIndex+i_maxLead; i++) {
				rProfiles << d_hillArr[i]/i_hillCount << ", ";
			}
			rProfiles << d_hillArr[i_maxIndex+i_maxLead]/i_hillCount << ")\n";


			rProfiles << "hillFallTime <- c(";
			for (int i = i_maxIndex-i_maxTrail; i < i_bamFiles-1; i++) {
				rProfiles << turnoverTimes[i] << ", ";
			}
			rProfiles << turnoverTimes[i_bamFiles-1] << ")\n";

			rProfiles << "hillFallDepth <- c(";
			for (int i = i_maxIndex-i_maxTrail; i < i_bamFiles-1; i++) {
				rProfiles << d_hillArr[i]/i_hillCount << ", ";
			}
			rProfiles << d_hillArr[i_bamFiles-1]/i_hillCount << ")\n";


			rProfiles << "hillTime <- c(";
			for (int i = 0; i < i_bamFiles-1; i++) {
				rProfiles << turnoverTimes[i] << ", ";
			}
			rProfiles << turnoverTimes[i_bamFiles-1] << ")\n";

			rProfiles << "hillDepth <- c(";
			for (int i = 0; i < i_bamFiles-1; i++) {
				rProfiles << d_hillArr[i]/i_hillCount << ", ";
			}
			rProfiles << d_hillArr[i_bamFiles-1]/i_hillCount << ")\n";
			rProfiles << "hill <- data.frame(hillTime,hillDepth)\n";
		}

		if (i_valley > 0) {
			//VALLEYS
			rProfiles << "valleyFallTime <- c(";
			for (int i = 0; i < i_minIndex+i_minLead; i++) {
				rProfiles << turnoverTimes[i] << ", ";
			}
			rProfiles << turnoverTimes[i_minIndex+i_minLead] << ")\n";

			rProfiles << "valleyFallDepth <- c(";
			for (int i = 0; i < i_minIndex+i_minLead; i++) {
				rProfiles << d_valleyArr[i]/i_valleyCount << ", ";
			}
			rProfiles << d_valleyArr[i_minIndex+i_minLead]/i_valleyCount << ")\n";


			rProfiles << "valleyRiseTime <- c(";
			for (int i = i_minIndex-i_minTrail; i < i_bamFiles-1; i++) {
				rProfiles << turnoverTimes[i] << ", ";
			}
			rProfiles << turnoverTimes[i_bamFiles-1] << ")\n";

			rProfiles << "valleyRiseDepth <- c(";
			for (int i = i_minIndex-i_minTrail; i < i_bamFiles-1; i++) {
				rProfiles << d_valleyArr[i]/i_valleyCount << ", ";
			}
			rProfiles << d_valleyArr[i_bamFiles-1]/i_valleyCount << ")\n";


			rProfiles << "valleyTime <- c(";
			for (int i = 0; i < i_bamFiles-1; i++) {
				rProfiles << turnoverTimes[i] << ", ";
			}
			rProfiles << turnoverTimes[i_bamFiles-1] << ")\n";

			rProfiles << "valleyDepth <- c(";
			for (int i = 0; i < i_bamFiles-1; i++) {
				rProfiles << d_valleyArr[i]/i_valleyCount << ", ";
			}
			rProfiles << d_valleyArr[i_bamFiles-1]/i_valleyCount << ")\n";
			rProfiles << "valley <- data.frame(valleyTime,valleyDepth)\n";
		}

		if (i_undefRise > 0) {
			// UNDEFINED RISES
			rProfiles << "undefRiseTime <- c(";
			for (int i = 0; i < i_bamFiles-1; i++) {
				rProfiles << turnoverTimes[i] << ", ";
			}
			rProfiles << turnoverTimes[i_bamFiles-1] << ")\n";

			rProfiles << "undefRiseDepth <- c(";
			for (int i = 0; i < i_bamFiles-1; i++) {
				rProfiles << d_undefRiseArr[i]/i_undefRiseCount << ", ";
			}
			rProfiles << d_undefRiseArr[i_bamFiles-1]/i_undefRiseCount << ")\n";
		}

		if (i_undefFall > 0) {
			// UNDEFINED FALLS
			rProfiles << "undefFallTime <- c(";
			for (int i = 0; i < i_bamFiles-1; i++) {
				rProfiles << turnoverTimes[i] << ", ";
			}
			rProfiles << turnoverTimes[i_bamFiles-1] << ")\n";

			rProfiles << "undefFallDepth <- c(";
			for (int i = 0; i < i_bamFiles-1; i++) {
				rProfiles << d_undefFallArr[i]/i_undefFallCount << ", ";
			}
			rProfiles << d_undefFallArr[i_bamFiles-1]/i_undefFallCount << ")\n";
		}
	}

	if (i_rise > 0) {
		//RISES
		rProfiles << "riseTime <- c(";
		for (int i = 0; i < i_bamFiles-1; i++) {
			rProfiles << turnoverTimes[i] << ", ";
		}
		rProfiles << turnoverTimes[i_bamFiles-1] << ")\n";

		rProfiles << "riseDepth <- c(";
		for (int i = 0; i < i_bamFiles-1; i++) {
			rProfiles << d_riseArr[i]/i_riseCount << ", ";
		}
		rProfiles << d_riseArr[i_bamFiles-1]/i_riseCount << ")\n";
	}
	if (i_fall > 0) {
		// FALLS
		rProfiles << "fallTime <- c(";
		for (int i = 0; i < i_bamFiles-1; i++) {
			rProfiles << turnoverTimes[i] << ", ";
		}
		rProfiles << turnoverTimes[i_bamFiles-1] << ")\n";

		rProfiles << "fallDepth <- c(";
		for (int i = 0; i < i_bamFiles-1; i++) {
			rProfiles << d_fallArr[i]/i_fallCount << ", ";
		}
		rProfiles << d_fallArr[i_bamFiles-1]/i_fallCount << ")\n";
	}

	if (i_rise > 0) {
		// RISES =======================================================================
		rProfiles << "Rise <- data.frame(riseTime,riseDepth)\n";
		rProfiles << "Rise.1 <- reshape2::melt(Rise,id.vars = \"riseTime\") # get numbers ready for use.\n";
		rProfiles << "Rise.fits <- expand.grid(conc2=seq(" << turnoverTimes[0] << ", " << turnoverTimes[i_bamFiles-1] << ", length=100))\n"; 

		// if else for lower asymptote below 0
		rProfiles << "success <- try(Rise.L.5 <- drm(riseDepth ~ riseTime, data = Rise, fct = L.5()),silent = TRUE)\n";
		rProfiles << "try(Rise.L.5f <- drm(riseDepth ~ riseTime, data = Rise, fct = L.5(fixed=c(NA,0,NA,NA,NA))),silent = TRUE)\n";

		rProfiles << "if (inherits(success,'try-error'))\n";
		rProfiles << " suppressWarnings(pm <- predict(Rise.L.5f, newdata=Rise.fits, interval=\"confidence\")) else\n";
		rProfiles << "    suppressWarnings(pm <- predict(Rise.L.5, newdata=Rise.fits, interval=\"confidence\")) \n";

		rProfiles << "Rise.fits$t <- pm[,1]\n";
		rProfiles << "Rise.fits$tmin <- pm[,2]\n";
		rProfiles << "Rise.fits$tmax <- pm[,3]\n";
		rProfiles << "riseProfile <- ggplot(Rise.1, aes(x = riseTime, y = riseDepth)) +\n";
		rProfiles << "	geom_point() +\n";
		rProfiles << "	#geom_ribbon(data=hillRise.fits, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2) +\n";
		rProfiles << "	geom_line(data=Rise.fits, aes(x=conc2, y=t), colour=\"chartreuse4\") +\n";
		rProfiles << "	labs(title = \"Average Rise Profile\") +\n";
		rProfiles << "	labs(x = \"Time\",size = 10) +\n";
		rProfiles << "	labs(y = \"Normalized Depth (reads)\",size = 10) +\n";
		for (int i = 0; i < i_bamFiles; i++) {
			rProfiles << "	geom_vline(xintercept = " << turnoverTimes[i] << ", colour=\"black\", alpha = 0.25, linetype = \"longdash\") +\n";
		}
		rProfiles << "	theme( \n";
		rProfiles << "		axis.text = element_text(size = 8, colour = \"black\"),\n";
		rProfiles << "		legend.key = element_rect(fill = \"white\"), \n";
		rProfiles << "		legend.background = element_rect(fill = \"white\"),\n";
		rProfiles << "		panel.grid.major = element_line(colour = \"white\"), \n";
		rProfiles << "		legend.position = \"bottom\", \n";
		rProfiles << "		panel.grid.minor = element_blank(),\n";
		rProfiles << "		panel.background = element_rect(fill = \"white\"),\n";
		rProfiles << "		axis.line.x = element_line(colour = \"black\"),\n";
		rProfiles << "		axis.line.y = element_line(colour = \"black\") \n";
		rProfiles << "		) \n";
	} else {
		rProfiles << "#VARIABLE T\n";
		rProfiles << "d = data.frame(x1=c(1), x2=c(2), y1=c(1), y2=c(2), t=c('No rises'))\n";
		rProfiles << "#VARIABLE NAME\n";
		rProfiles << "riseProfile <- ggplot() + \n";
		rProfiles << "scale_y_continuous(expand = c(0,0), limits = c(0, 2.5)) + \n";
		rProfiles << "scale_x_continuous(expand = c(0,0), limits = c(-2, 6)) + \n";
		rProfiles << "geom_text(data=d, aes(x=1.5, y=1.5, label=t), size=6) +\n";
		rProfiles << "	theme( \n";
		rProfiles << "		axis.text = element_text(size = 8, colour = \"white\"),\n";
		rProfiles << "		legend.key = element_rect(fill = \"white\"), \n";
		rProfiles << "		legend.background = element_rect(fill = \"white\"),\n";
		rProfiles << "		panel.grid.major = element_line(colour = \"white\"), \n";
		rProfiles << "		legend.position = \"\", \n";
		rProfiles << "		panel.grid.minor = element_blank(),\n";
		rProfiles << "		panel.background = element_rect(fill = \"white\"),\n";
		rProfiles << "		axis.line.x = element_line(colour = \"white\"),\n";
		rProfiles << "		axis.line.y = element_line(colour = \"white\"),\n";
		rProfiles << "		axis.title.x=element_blank(),\n";
		rProfiles << "		axis.text.x=element_blank(),\n";
		rProfiles << "		axis.ticks.x=element_blank(),\n";
		rProfiles << "		axis.title.y=element_blank(),\n";
		rProfiles << "		axis.text.y=element_blank(),\n";
		rProfiles << "		axis.ticks.y=element_blank()\n";
		rProfiles << "		) \n";
	}
	if (i_fall > 0) {

		// FALLS ===========================================================================
		rProfiles << "fall <- data.frame(fallTime,fallDepth)\n";
		rProfiles << "fall.1 <- reshape2::melt(fall,id.vars = \"fallTime\") # get numbers ready for use.\n";
		rProfiles << "fall.fits <- expand.grid(conc2=seq(" << turnoverTimes[0] << ", " << turnoverTimes[i_bamFiles-1] << ", length=100))\n"; 

		// if else for lower asymptote below 0
		rProfiles << "success <- try(fall.L.5 <- drm(fallDepth ~ fallTime, data = fall, fct = L.5()),silent = TRUE)\n";
		rProfiles << "try(fall.L.5f <- drm(fallDepth ~ fallTime, data = fall, fct = L.5(fixed=c(NA,0,NA,NA,NA))),silent = TRUE)\n";

		rProfiles << "if (inherits(success,'try-error'))\n";
		rProfiles << " suppressWarnings(pm <- predict(fall.L.5f, newdata=fall.fits, interval=\"confidence\")) else\n";
		rProfiles << "    suppressWarnings(pm <- predict(fall.L.5, newdata=fall.fits, interval=\"confidence\")) \n";

		rProfiles << "fall.fits$t <- pm[,1]\n";
		rProfiles << "fall.fits$tmin <- pm[,2]\n";
		rProfiles << "fall.fits$tmax <- pm[,3]\n";
		rProfiles << "fallProfile <- ggplot(fall.1, aes(x = fallTime, y = fallDepth)) +\n";
		rProfiles << "	geom_point() +\n";
		rProfiles << "	#geom_ribbon(data=hillfall.fits, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2) +\n";
		rProfiles << "	geom_line(data=fall.fits, aes(x=conc2, y=t), colour=\"chartreuse4\") +\n";
		rProfiles << "	labs(title = \"Average Fall Profile\") +\n";
		rProfiles << "	labs(x = \"Time\",size = 10) +\n";
		rProfiles << "	labs(y = \"Normalized Depth (reads)\",size = 10) +\n";
		for (int i = 0; i < i_bamFiles; i++) {
			rProfiles << "	geom_vline(xintercept = " << turnoverTimes[i] << ", colour=\"black\", alpha = 0.25, linetype = \"longdash\") +\n";
		}
		rProfiles << "	theme( \n";
		rProfiles << "		axis.text = element_text(size = 8, colour = \"black\"),\n";
		rProfiles << "		legend.key = element_rect(fill = \"white\"), \n";
		rProfiles << "		legend.background = element_rect(fill = \"white\"),\n";
		rProfiles << "		panel.grid.major = element_line(colour = \"white\"), \n";
		rProfiles << "		legend.position = \"bottom\", \n";
		rProfiles << "		panel.grid.minor = element_blank(),\n";
		rProfiles << "		panel.background = element_rect(fill = \"white\"),\n";
		rProfiles << "		axis.line.x = element_line(colour = \"black\"),\n";
		rProfiles << "		axis.line.y = element_line(colour = \"black\") \n";
		rProfiles << "		) \n";
	} else {
		rProfiles << "#VARIABLE T\n";
		rProfiles << "d = data.frame(x1=c(1), x2=c(2), y1=c(1), y2=c(2), t=c('No falls'))\n";
		rProfiles << "#VARIABLE NAME\n";
		rProfiles << "fallProfile <- ggplot() + \n";
		rProfiles << "scale_y_continuous(expand = c(0,0), limits = c(0, 2.5)) + \n";
		rProfiles << "scale_x_continuous(expand = c(0,0), limits = c(-2, 6)) + \n";
		rProfiles << "geom_text(data=d, aes(x=1.5, y=1.5, label=t), size=6) +\n";
		rProfiles << "	theme( \n";
		rProfiles << "		axis.text = element_text(size = 8, colour = \"white\"),\n";
		rProfiles << "		legend.key = element_rect(fill = \"white\"), \n";
		rProfiles << "		legend.background = element_rect(fill = \"white\"),\n";
		rProfiles << "		panel.grid.major = element_line(colour = \"white\"), \n";
		rProfiles << "		legend.position = \"\", \n";
		rProfiles << "		panel.grid.minor = element_blank(),\n";
		rProfiles << "		panel.background = element_rect(fill = \"white\"),\n";
		rProfiles << "		axis.line.x = element_line(colour = \"white\"),\n";
		rProfiles << "		axis.line.y = element_line(colour = \"white\"),\n";
		rProfiles << "		axis.title.x=element_blank(),\n";
		rProfiles << "		axis.text.x=element_blank(),\n";
		rProfiles << "		axis.ticks.x=element_blank(),\n";
		rProfiles << "		axis.title.y=element_blank(),\n";
		rProfiles << "		axis.text.y=element_blank(),\n";
		rProfiles << "		axis.ticks.y=element_blank()\n";
		rProfiles << "		) \n";
	}



	if (b_model) {
		if (i_hill > 0) {

			// HILL ===============================================================
			rProfiles << "hillRise <- data.frame(hillRiseTime,hillRiseDepth)\n";
			rProfiles << "hillRise.1 <- reshape2::melt(hillRise,id.vars = \"hillRiseTime\") # get numbers ready for use.\n";
			rProfiles << "hillRise.fits <- expand.grid(conc=seq(" << turnoverTimes[0] << ", " << turnoverTimes[i_bamFiles-1] << ", length=100))\n"; 

			// if else for lower asymptote below 0
			rProfiles << "success <- try(hillRise.L.5 <- drm(hillRiseDepth ~ hillRiseTime, data = hillRise, fct = L.5()),silent = TRUE)\n";
			rProfiles << "try(hillRise.L.5f <- drm(hillRiseDepth ~ hillRiseTime, data = hillRise, fct = L.5(fixed=c(NA,0,NA,NA,NA))),silent = TRUE)\n";

			rProfiles << "if (inherits(success,'try-error'))\n";
			rProfiles << " suppressWarnings(pm <- predict(hillRise.L.5f, newdata=hillRise.fits, interval=\"confidence\")) else\n";
			rProfiles << "    suppressWarnings(pm <- predict(hillRise.L.5, newdata=hillRise.fits, interval=\"confidence\")) \n";

			rProfiles << "hillRise.fits$t <- pm[,1]\n";
			rProfiles << "hillRise.fits$tmin <- pm[,2]\n";
			rProfiles << "hillRise.fits$tmax <- pm[,3]\n";


			rProfiles << "hillFall <- data.frame(hillFallTime,hillFallDepth)\n";
			rProfiles << "hillFall.1 <- reshape2::melt(hillFall,id.vars = \"hillFallTime\") # get numbers ready for use.\n";
			rProfiles << "hillFall.fits <- expand.grid(conc2=seq(" << turnoverTimes[0] << ", " << turnoverTimes[i_bamFiles-1] << ", length=100))\n"; 

			// if else for lower asymptote below 0
			rProfiles << "success <- try(hillFall.L.5 <- drm(hillFallDepth ~ hillFallTime, data = hillFall, fct = L.5()),silent = TRUE)\n";
			rProfiles << "try(hillFall.L.5f <- drm(hillFallDepth ~ hillFallTime, data = hillFall, fct = L.5(fixed=c(NA,0,NA,NA,NA))),silent = TRUE)\n";

			rProfiles << "if (inherits(success,'try-error'))\n";
			rProfiles << " suppressWarnings(pm <- predict(hillFall.L.5f, newdata=hillFall.fits, interval=\"confidence\")) else\n";
			rProfiles << "    suppressWarnings(pm <- predict(hillFall.L.5, newdata=hillFall.fits, interval=\"confidence\")) \n";

			rProfiles << "hillFall.fits$t <- pm[,1]\n";
			rProfiles << "hillFall.fits$tmin <- pm[,2]\n";
			rProfiles << "hillFall.fits$tmax <- pm[,3]\n";

			rProfiles << "hillProfile <- ggplot(hill, aes(x = hillTime, y = hillDepth)) +\n";
			rProfiles << "	geom_point() +\n";
			rProfiles << "	geom_line(data=hillRise.fits, aes(x=conc, y=t), colour=\"navy\") +\n";
			rProfiles << "	geom_line(data=hillFall.fits, aes(x=conc2, y=t), colour=\"chartreuse4\") +\n";
			rProfiles << "	labs(title = \"Average Hill Profile\") +\n";
			rProfiles << "	labs(x = \"Time\",size = 10) +\n";
			rProfiles << "	labs(y = \"Normalized Depth (reads)\",size = 10) +\n";
			for (int i = 0; i < i_bamFiles; i++) {
				rProfiles << "	geom_vline(xintercept = " << turnoverTimes[i] << ", colour=\"black\", alpha = 0.25, linetype = \"longdash\") +\n";
			}
			rProfiles << "	theme( \n";
			rProfiles << "		axis.text = element_text(size = 8, colour = \"black\"),\n";
			rProfiles << "		legend.key = element_rect(fill = \"white\"), \n";
			rProfiles << "		legend.background = element_rect(fill = \"white\"),\n";
			rProfiles << "		panel.grid.major = element_line(colour = \"white\"), \n";
			rProfiles << "		legend.position = \"bottom\", \n";
			rProfiles << "		panel.grid.minor = element_blank(),\n";
			rProfiles << "		panel.background = element_rect(fill = \"white\"),\n";
			rProfiles << "		axis.line.x = element_line(colour = \"black\"),\n";
			rProfiles << "		axis.line.y = element_line(colour = \"black\") \n";
			rProfiles << "		) \n";
		} else {
			rProfiles << "#VARIABLE T\n";
			rProfiles << "d = data.frame(x1=c(1), x2=c(2), y1=c(1), y2=c(2), t=c('No hills'))\n";
			rProfiles << "#VARIABLE NAME\n";
			rProfiles << "hillProfile <- ggplot() + \n";
			rProfiles << "scale_y_continuous(expand = c(0,0), limits = c(0, 2.5)) + \n";
			rProfiles << "scale_x_continuous(expand = c(0,0), limits = c(-2, 6)) + \n";
			rProfiles << "geom_text(data=d, aes(x=1.5, y=1.5, label=t), size=6) +\n";
			rProfiles << "	theme( \n";
			rProfiles << "		axis.text = element_text(size = 8, colour = \"white\"),\n";
			rProfiles << "		legend.key = element_rect(fill = \"white\"), \n";
			rProfiles << "		legend.background = element_rect(fill = \"white\"),\n";
			rProfiles << "		panel.grid.major = element_line(colour = \"white\"), \n";
			rProfiles << "		legend.position = \"\", \n";
			rProfiles << "		panel.grid.minor = element_blank(),\n";
			rProfiles << "		panel.background = element_rect(fill = \"white\"),\n";
			rProfiles << "		axis.line.x = element_line(colour = \"white\"),\n";
			rProfiles << "		axis.line.y = element_line(colour = \"white\"),\n";
			rProfiles << "		axis.title.x=element_blank(),\n";
			rProfiles << "		axis.text.x=element_blank(),\n";
			rProfiles << "		axis.ticks.x=element_blank(),\n";
			rProfiles << "		axis.title.y=element_blank(),\n";
			rProfiles << "		axis.text.y=element_blank(),\n";
			rProfiles << "		axis.ticks.y=element_blank()\n";
			rProfiles << "		) \n";
		}
		if (i_valley > 0) {
			// VALLEYS =============================================================

				rProfiles << "valleyRise <- data.frame(valleyRiseTime,valleyRiseDepth)\n";
			rProfiles << "valleyRise.1 <- reshape2::melt(valleyRise,id.vars = \"valleyRiseTime\") # get numbers ready for use.\n";
			rProfiles << "valleyRise.fits <- expand.grid(conc=seq(" << turnoverTimes[0] << ", " << turnoverTimes[i_bamFiles-1] << ", length=100))\n"; 

			// if else for lower asymptote below 0
			rProfiles << "success <- try(valleyRise.L.5 <- drm(valleyRiseDepth ~ valleyRiseTime, data = valleyRise, fct = L.5()),silent = TRUE)\n";
			rProfiles << "try(valleyRise.L.5f <- drm(valleyRiseDepth ~ valleyRiseTime, data = valleyRise, fct = L.5(fixed=c(NA,0,NA,NA,NA))),silent = TRUE)\n";

			rProfiles << "if (inherits(success,'try-error'))\n";
			rProfiles << " suppressWarnings(pm <- predict(valleyRise.L.5f, newdata=valleyRise.fits, interval=\"confidence\")) else\n";
			rProfiles << "    suppressWarnings(pm <- predict(valleyRise.L.5, newdata=valleyRise.fits, interval=\"confidence\")) \n";

			rProfiles << "valleyRise.fits$t <- pm[,1]\n";
			rProfiles << "valleyRise.fits$tmin <- pm[,2]\n";
			rProfiles << "valleyRise.fits$tmax <- pm[,3]\n";


			rProfiles << "valleyFall <- data.frame(valleyFallTime,valleyFallDepth)\n";
			rProfiles << "valleyFall.1 <- reshape2::melt(valleyFall,id.vars = \"valleyFallTime\") # get numbers ready for use.\n";
			rProfiles << "valleyFall.fits <- expand.grid(conc2=seq(" << turnoverTimes[0] << ", " << turnoverTimes[i_bamFiles-1] << ", length=100))\n"; 

			// if else for lower asymptote below 0
			rProfiles << "success <- try(valleyFall.L.5 <- drm(valleyFallDepth ~ valleyFallTime, data = valleyFall, fct = L.5()),silent = TRUE)\n";
			rProfiles << "try(valleyFall.L.5f <- drm(valleyFallDepth ~ valleyFallTime, data = valleyFall, fct = L.5(fixed=c(NA,0,NA,NA,NA))),silent = TRUE)\n";

			rProfiles << "if (inherits(success,'try-error'))\n";
			rProfiles << " suppressWarnings(pm <- predict(valleyFall.L.5f, newdata=valleyFall.fits, interval=\"confidence\")) else\n";
			rProfiles << "    suppressWarnings(pm <- predict(valleyFall.L.5, newdata=valleyFall.fits, interval=\"confidence\")) \n";

			rProfiles << "valleyFall.fits$t <- pm[,1]\n";
			rProfiles << "valleyFall.fits$tmin <- pm[,2]\n";
			rProfiles << "valleyFall.fits$tmax <- pm[,3]\n";

			rProfiles << "valleyProfile <- ggplot(valley, aes(x = valleyTime, y = valleyDepth)) +\n";
			rProfiles << "	geom_point() +\n";
			rProfiles << "	geom_line(data=valleyRise.fits, aes(x=conc, y=t), colour=\"navy\") +\n";
			rProfiles << "	geom_line(data=valleyFall.fits, aes(x=conc2, y=t), colour=\"chartreuse4\") +\n";
			rProfiles << "	labs(title = \"Average Valley Profile\") +\n";
			rProfiles << "	labs(x = \"Time\",size = 10) +\n";
			rProfiles << "	labs(y = \"Normalized Depth (reads)\",size = 10) +\n";
			for (int i = 0; i < i_bamFiles; i++) {
				rProfiles << "	geom_vline(xintercept = " << turnoverTimes[i] << ", colour=\"black\", alpha = 0.25, linetype = \"longdash\") +\n";
			}
			rProfiles << "	theme( \n";
			rProfiles << "		axis.text = element_text(size = 8, colour = \"black\"),\n";
			rProfiles << "		legend.key = element_rect(fill = \"white\"), \n";
			rProfiles << "		legend.background = element_rect(fill = \"white\"),\n";
			rProfiles << "		panel.grid.major = element_line(colour = \"white\"), \n";
			rProfiles << "		legend.position = \"bottom\", \n";
			rProfiles << "		panel.grid.minor = element_blank(),\n";
			rProfiles << "		panel.background = element_rect(fill = \"white\"),\n";
			rProfiles << "		axis.line.x = element_line(colour = \"black\"),\n";
			rProfiles << "		axis.line.y = element_line(colour = \"black\") \n";
			rProfiles << "		) \n";
		} else {
			rProfiles << "#VARIABLE T\n";
			rProfiles << "d = data.frame(x1=c(1), x2=c(2), y1=c(1), y2=c(2), t=c('No vallies'))\n";
			rProfiles << "#VARIABLE NAME\n";
			rProfiles << "valleyProfile <- ggplot() + \n";
			rProfiles << "scale_y_continuous(expand = c(0,0), limits = c(0, 2.5)) + \n";
			rProfiles << "scale_x_continuous(expand = c(0,0), limits = c(-2, 6)) + \n";
			rProfiles << "geom_text(data=d, aes(x=1.5, y=1.5, label=t), size=6) +\n";
			rProfiles << "	theme( \n";
			rProfiles << "		axis.text = element_text(size = 8, colour = \"white\"),\n";
			rProfiles << "		legend.key = element_rect(fill = \"white\"), \n";
			rProfiles << "		legend.background = element_rect(fill = \"white\"),\n";
			rProfiles << "		panel.grid.major = element_line(colour = \"white\"), \n";
			rProfiles << "		legend.position = \"\", \n";
			rProfiles << "		panel.grid.minor = element_blank(),\n";
			rProfiles << "		panel.background = element_rect(fill = \"white\"),\n";
			rProfiles << "		axis.line.x = element_line(colour = \"white\"),\n";
			rProfiles << "		axis.line.y = element_line(colour = \"white\"),\n";
			rProfiles << "		axis.title.x=element_blank(),\n";
			rProfiles << "		axis.text.x=element_blank(),\n";
			rProfiles << "		axis.ticks.x=element_blank(),\n";
			rProfiles << "		axis.title.y=element_blank(),\n";
			rProfiles << "		axis.text.y=element_blank(),\n";
			rProfiles << "		axis.ticks.y=element_blank()\n";
			rProfiles << "		) \n";
		}
		if (i_undefRise > 0) {
			// UNDEFINED RISES =======================================================
			rProfiles << "undefRise <- data.frame(undefRiseTime,undefRiseDepth)\n";
			rProfiles << "undefRise.1 <- reshape2::melt(undefRise,id.vars = \"undefRiseTime\") # get numbers ready for use.\n";
			rProfiles << "undefRise.fits <- expand.grid(conc2=seq(" << turnoverTimes[0] << ", " << turnoverTimes[i_bamFiles-1] << ", length=100))\n"; 

			// if else for lower asymptote below 0
			rProfiles << "success <- try(undefRise.L.5 <- drm(undefRiseDepth ~ undefRiseTime, data = undefRise, fct = L.5()),silent = TRUE)\n";
			rProfiles << "try(undefRise.L.5f <- drm(undefRiseDepth ~ undefRiseTime, data = undefRise, fct = L.5(fixed=c(NA,0,NA,NA,NA))),silent = TRUE)\n";

			rProfiles << "if (inherits(success,'try-error'))\n";
			rProfiles << " suppressWarnings(pm <- predict(undefRise.L.5f, newdata=undefRise.fits, interval=\"confidence\")) else\n";
			rProfiles << "    suppressWarnings(pm <- predict(undefRise.L.5, newdata=undefRise.fits, interval=\"confidence\")) \n";

			rProfiles << "undefRise.fits$t <- pm[,1]\n";
			rProfiles << "undefRise.fits$tmin <- pm[,2]\n";
			rProfiles << "undefRise.fits$tmax <- pm[,3]\n";
			rProfiles << "undefRiseProfile <- ggplot(undefRise.1, aes(x = undefRiseTime, y = undefRiseDepth)) +\n";
			rProfiles << "	geom_point() +\n";
			rProfiles << "	#geom_ribbon(data=hillRise.fits, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2) +\n";
			rProfiles << "	geom_line(data=undefRise.fits, aes(x=conc2, y=t), colour=\"chartreuse4\") +\n";
			rProfiles << "	labs(title = \"Average Undefined Rise Profile\") +\n";
			rProfiles << "	labs(x = \"Time\",size = 10) +\n";
			rProfiles << "	labs(y = \"Normalized Depth (reads)\",size = 10) +\n";
			for (int i = 0; i < i_bamFiles; i++) {
				rProfiles << "	geom_vline(xintercept = " << turnoverTimes[i] << ", colour=\"black\", alpha = 0.25, linetype = \"longdash\") +\n";
			}
			rProfiles << "	theme( \n";
			rProfiles << "		axis.text = element_text(size = 8, colour = \"black\"),\n";
			rProfiles << "		legend.key = element_rect(fill = \"white\"), \n";
			rProfiles << "		legend.background = element_rect(fill = \"white\"),\n";
			rProfiles << "		panel.grid.major = element_line(colour = \"white\"), \n";
			rProfiles << "		legend.position = \"bottom\", \n";
			rProfiles << "		panel.grid.minor = element_blank(),\n";
			rProfiles << "		panel.background = element_rect(fill = \"white\"),\n";
			rProfiles << "		axis.line.x = element_line(colour = \"black\"),\n";
			rProfiles << "		axis.line.y = element_line(colour = \"black\") \n";
			rProfiles << "		) \n";
		} else {
			rProfiles << "#VARIABLE T\n";
			rProfiles << "d = data.frame(x1=c(1), x2=c(2), y1=c(1), y2=c(2), t=c('No undefined rises'))\n";
			rProfiles << "#VARIABLE NAME\n";
			rProfiles << "undefRiseProfile <- ggplot() + \n";
			rProfiles << "scale_y_continuous(expand = c(0,0), limits = c(0, 2.5)) + \n";
			rProfiles << "scale_x_continuous(expand = c(0,0), limits = c(-2, 6)) + \n";
			rProfiles << "geom_text(data=d, aes(x=1.5, y=1.5, label=t), size=6) +\n";
			rProfiles << "	theme( \n";
			rProfiles << "		axis.text = element_text(size = 8, colour = \"white\"),\n";
			rProfiles << "		legend.key = element_rect(fill = \"white\"), \n";
			rProfiles << "		legend.background = element_rect(fill = \"white\"),\n";
			rProfiles << "		panel.grid.major = element_line(colour = \"white\"), \n";
			rProfiles << "		legend.position = \"\", \n";
			rProfiles << "		panel.grid.minor = element_blank(),\n";
			rProfiles << "		panel.background = element_rect(fill = \"white\"),\n";
			rProfiles << "		axis.line.x = element_line(colour = \"white\"),\n";
			rProfiles << "		axis.line.y = element_line(colour = \"white\"),\n";
			rProfiles << "		axis.title.x=element_blank(),\n";
			rProfiles << "		axis.text.x=element_blank(),\n";
			rProfiles << "		axis.ticks.x=element_blank(),\n";
			rProfiles << "		axis.title.y=element_blank(),\n";
			rProfiles << "		axis.text.y=element_blank(),\n";
			rProfiles << "		axis.ticks.y=element_blank()\n";
			rProfiles << "		) \n";
		}
		if (i_undefFall > 0) {
			// UNDEFINED FALLS ==============================================================
			rProfiles << "undefFall <- data.frame(undefFallTime,undefFallDepth)\n";
			rProfiles << "undefFall.1 <- reshape2::melt(undefFall,id.vars = \"undefFallTime\") # get numbers ready for use.\n";
			rProfiles << "undefFall.fits <- expand.grid(conc2=seq(" << turnoverTimes[0] << ", " << turnoverTimes[i_bamFiles-1] << ", length=100))\n"; 

			// if else for lower asymptote below 0
			rProfiles << "success <- try(undefFall.L.5 <- drm(undefFallDepth ~ undefFallTime, data = undefFall, fct = L.5()),silent = TRUE)\n";
			rProfiles << "try(undefFall.L.5f <- drm(undefFallDepth ~ undefFallTime, data = undefFall, fct = L.5(fixed=c(NA,0,NA,NA,NA))),silent = TRUE)\n";

			rProfiles << "if (inherits(success,'try-error'))\n";
			rProfiles << " suppressWarnings(pm <- predict(undefFall.L.5f, newdata=undefFall.fits, interval=\"confidence\")) else\n";
			rProfiles << "    suppressWarnings(pm <- predict(undefFall.L.5, newdata=undefFall.fits, interval=\"confidence\")) \n";

			rProfiles << "undefFall.fits$t <- pm[,1]\n";
			rProfiles << "undefFall.fits$tmin <- pm[,2]\n";
			rProfiles << "undefFall.fits$tmax <- pm[,3]\n";
			rProfiles << "undefFallProfile <- ggplot(undefFall.1, aes(x = undefFallTime, y = undefFallDepth)) +\n";
			rProfiles << "	geom_point() +\n";
			rProfiles << "	#geom_ribbon(data=hillundefFall.fits, aes(x=conc, y=p, ymin=pmin, ymax=pmax), alpha=0.2) +\n";
			rProfiles << "	geom_line(data=undefFall.fits, aes(x=conc2, y=t), colour=\"chartreuse4\") +\n";
			rProfiles << "	labs(title = \"Average Undefined Fall Profile\") +\n";
			rProfiles << "	labs(x = \"Time\",size = 10) +\n";
			rProfiles << "	labs(y = \"Normalized Depth (reads)\",size = 10) +\n";
			for (int i = 0; i < i_bamFiles; i++) {
				rProfiles << "	geom_vline(xintercept = " << turnoverTimes[i] << ", colour=\"black\", alpha = 0.25, linetype = \"longdash\") +\n";
			}
			rProfiles << "	theme( \n";
			rProfiles << "		axis.text = element_text(size = 8, colour = \"black\"),\n";
			rProfiles << "		legend.key = element_rect(fill = \"white\"), \n";
			rProfiles << "		legend.background = element_rect(fill = \"white\"),\n";
			rProfiles << "		panel.grid.major = element_line(colour = \"white\"), \n";
			rProfiles << "		legend.position = \"bottom\", \n";
			rProfiles << "		panel.grid.minor = element_blank(),\n";
			rProfiles << "		panel.background = element_rect(fill = \"white\"),\n";
			rProfiles << "		axis.line.x = element_line(colour = \"black\"),\n";
			rProfiles << "		axis.line.y = element_line(colour = \"black\") \n";
			rProfiles << "		) \n";
		} else {
			rProfiles << "#VARIABLE T\n";
			rProfiles << "d = data.frame(x1=c(1), x2=c(2), y1=c(1), y2=c(2), t=c('No undefined falls'))\n";
			rProfiles << "#VARIABLE NAME\n";
			rProfiles << "undefFallProfile <- ggplot() + \n";
			rProfiles << "scale_y_continuous(expand = c(0,0), limits = c(0, 2.5)) + \n";
			rProfiles << "scale_x_continuous(expand = c(0,0), limits = c(-2, 6)) + \n";
			rProfiles << "geom_text(data=d, aes(x=1.5, y=1.5, label=t), size=6) +\n";
			rProfiles << "	theme( \n";
			rProfiles << "		axis.text = element_text(size = 8, colour = \"white\"),\n";
			rProfiles << "		legend.key = element_rect(fill = \"white\"), \n";
			rProfiles << "		legend.background = element_rect(fill = \"white\"),\n";
			rProfiles << "		panel.grid.major = element_line(colour = \"white\"), \n";
			rProfiles << "		legend.position = \"\", \n";
			rProfiles << "		panel.grid.minor = element_blank(),\n";
			rProfiles << "		panel.background = element_rect(fill = \"white\"),\n";
			rProfiles << "		axis.line.x = element_line(colour = \"white\"),\n";
			rProfiles << "		axis.line.y = element_line(colour = \"white\"),\n";
			rProfiles << "		axis.title.x=element_blank(),\n";
			rProfiles << "		axis.text.x=element_blank(),\n";
			rProfiles << "		axis.ticks.x=element_blank(),\n";
			rProfiles << "		axis.title.y=element_blank(),\n";
			rProfiles << "		axis.text.y=element_blank(),\n";
			rProfiles << "		axis.ticks.y=element_blank()\n";
			rProfiles << "		) \n";
		}
		
	}



	rProfiles.close();
}












