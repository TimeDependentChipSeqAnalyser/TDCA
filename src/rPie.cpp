#include <iostream>
#include <fstream>
#include <string>
#include <vector>
using namespace std;

void rPie(std::string s_rPltsName, int i_peakNumber, int i_undef, int i_elim, int i_valley, int i_hill, int i_fall, int i_rise, bool b_model)
{



	// write rthon script
	ofstream rPiePlt;
	rPiePlt.open(s_rPltsName, std::ios::app);

	double d_undef = double(i_undef)/double(i_peakNumber)*100;
	double d_elim = double(i_elim)/double(i_peakNumber)*100;
	double d_hill = double(i_hill)/double(i_peakNumber)*100;
	double d_rise = double(i_rise)/double(i_peakNumber)*100;
	double d_valley = double(i_valley)/double(i_peakNumber)*100;
	double d_fall = double(i_fall)/double(i_peakNumber)*100;

	if (b_model) {
		rPiePlt << "values <- c(" << d_valley << "," << d_undef << "," << d_rise << "," << d_hill << "," << d_fall << "," << d_elim << ")\n";
		rPiePlt << "id <- c(\"valley\", \"undefined\", \"rise\", \"hill\", \"fall\", \"eliminated\")\n";
	} else {
		rPiePlt << "values <- c(" << d_rise << "," << d_fall << ")\n";
		rPiePlt << "id <- c(\"rise\", \"fall\")\n";
	}

	rPiePlt << "dfvalue <- data.frame(values, id)\n";

	rPiePlt << "bp<- ggplot(dfvalue, aes(x=\"\", y=values, fill=id))+\n";
	rPiePlt << "  geom_bar(width = 1, stat = \"identity\", color='black')\n";
	rPiePlt << "pie <- bp + coord_polar(\"y\", start=0)\n";

	rPiePlt << "# coordinate midpoints\n";
	rPiePlt << "y.breaks <- cumsum(dfvalue$values) - dfvalue$values/2\n";


	rPiePlt << "blank_theme <- theme_minimal()+\n";
	rPiePlt << "	theme(\n";
	rPiePlt << "		axis.title.x = element_blank(),\n";
	rPiePlt << "		axis.title.y = element_blank(),\n";
	rPiePlt << "		panel.border = element_blank(),\n";
	rPiePlt << "		axis.text.x = element_text(size=6),\n";
	rPiePlt << "		panel.grid=element_blank(),\n";
	rPiePlt << "		axis.ticks = element_blank(),\n";
	rPiePlt << "		plot.title=element_text(size=10),\n";		
	rPiePlt << "		legend.key.size = unit(0.25, \"cm\"),\n";			
	rPiePlt << "		legend.text = element_text(size = 7)\n";
	rPiePlt << "	)\n";

	rPiePlt << "lbls <- paste(round(values,1),\"%\",sep=\"\")\n";
	if (b_model) {
		rPiePlt << "pie <- pie + scale_fill_manual(values=c(\"firebrick2\", \"seashell2\", \"khaki3\", \"khaki2\", \"darkorange\", \"seashell3\")) + blank_theme +\n";
	} else {
		rPiePlt << "pie <- pie + scale_fill_manual(values=c(\"khaki2\", \"seashell2\")) + blank_theme +\n";
	}
	rPiePlt << "	theme(axis.text.x=element_text(color='black')) +\n";
	rPiePlt << "	labs(fill=\"\") +\n";
	rPiePlt << "	guides(fill=guide_legend(nrow=2,byrow=TRUE)) +\n";
	rPiePlt << "	scale_y_continuous(breaks=y.breaks, labels=lbls) +\n";
	rPiePlt << "	ggtitle(\"Peak Characteristics\")\n";


	rPiePlt.close();

}
