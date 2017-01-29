#include <iostream>
#include <fstream>
#include <string>

// creates scatterplot boxplot of log2(rise/fall) of upper aymptotes in hills and valleys.
void rBlank(std::string s_rPltsName, std::string s_direction, std::string s_plotName) 
{

	// write R script
	std::ofstream rBlank;
	rBlank.open (s_rPltsName, std::ios::app);

	rBlank << "#VARIABLE T\n";
	rBlank << "d = data.frame(x1=c(1), x2=c(2), y1=c(1), y2=c(2), t=c('No " << s_direction << " data'))\n";
	rBlank << "#VARIABLE NAME\n";
	rBlank <<  s_plotName << " <- ggplot() + \n";
	rBlank << "scale_y_continuous(expand = c(0,0), limits = c(0, 2.5)) + \n";
	rBlank << "scale_x_continuous(expand = c(0,0), limits = c(-2, 6)) + \n";
	rBlank << "geom_text(data=d, aes(x=1.5, y=1.5, label=t), size=16) +\n";
	rBlank << "	theme( \n";
	rBlank << "		axis.text = element_text(size = 8, colour = \"white\"),\n";
	rBlank << "		legend.key = element_rect(fill = \"white\"), \n";
	rBlank << "		legend.background = element_rect(fill = \"white\"),\n";
	rBlank << "		panel.grid.major = element_line(colour = \"white\"), \n";
	rBlank << "		legend.position = \"\", \n";
	rBlank << "		panel.grid.minor = element_blank(),\n";
	rBlank << "		panel.background = element_rect(fill = \"white\"),\n";
	rBlank << "		axis.line.x = element_line(colour = \"white\"),\n";
	rBlank << "		axis.line.y = element_line(colour = \"white\"),\n";
	rBlank << "		axis.title.x=element_blank(),\n";
	rBlank << "		axis.text.x=element_blank(),\n";
	rBlank << "		axis.ticks.x=element_blank(),\n";
	rBlank << "		axis.title.y=element_blank(),\n";
	rBlank << "		axis.text.y=element_blank(),\n";
	rBlank << "		axis.ticks.y=element_blank()\n";
	rBlank << "		) \n";

	rBlank.close();
}


