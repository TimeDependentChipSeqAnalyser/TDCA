#include <iostream>
#include <sstream>
#include <string>

int welcome()
{
	std::cout << "Welcome to TDCA. For user manual and detailed usage see: https://github.com/TimeDependentChipSeqAnalyser/TDCA" << std::endl;
	std::cout << "Usage: tdca <-bed peaks.bed> <-bam folder/> [options]" << std::endl;
	std::cout << "Example: tdca -bed ChIP-seq.peaks.bed -bam bamFolder/ -i bamInputFolder/ -g mm9 -name exp-name" <<
		std::endl;
	std::cout << "Options:" << std::endl;
	std::cout << "    -v			Display program version and exit program." << std::endl;
	std::cout << "    -h			Display help page and exit program." << std::endl;
	std::cout << "    -bam 		User specified folder containing sorted bam turnover files " <<
							"including index files." << std::endl;
	std::cout << "    -i 			User specified folder containing sorted input bam turnover " << 
							"files including index files." << std::endl;
	std::cout << "    -bed 		User specified bed file containing loci of interest." << std::endl;
	std::cout << "    -g 			Genome name (hg19, hg38, mm9, mm10, dm3, dm6, ce11, sacCer3)." << std::endl;
	std::cout << "    -3d 		User specified gene file containing RefSeq gene names." << std::endl;
	std::cout << "    -s 			Saturation threshold consideration for data modelling (allowable range from 0.5-0.95). Default = 0.85." << std::endl;
	std::cout << "    -t 			Time point threshold consideration for data modelling (allowable range from 0-2). Default = 1." << std::endl;
	std::cout << "    -name	 	User specified name for output files." << std::endl;
	std::cout << "    -model	 	Data modelled based on prediction. Default is no data modelling." << std::endl;
	std::cout << "    -dm	 		Data matrix used to normalize user defined input files." << std::endl;
	std::cout << "    -p	 		Explicitly state number of processors to use." << std::endl;
	std::cout << "    -L4	 		Data modelled to 4 parameter logistic curve (no asymetry value)." << std::endl;
	return 0;
}
