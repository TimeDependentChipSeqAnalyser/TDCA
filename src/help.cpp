#include <iostream>
#include <sstream>
#include <string>

int help()
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
	std::cout << "    -L5	 		Data modelled to 5 parameter sigmoidal curve (asymmetry factor applies)." << std::endl;
	std::cout << "    -poisson	 	Data modelled to 3 parameter sigmoidal curve assuming Poisson distributed coverage." << std::endl;
	std::cout << "    -dm	 		Data matrix used to normalize user defined input files." << std::endl;
	std::cout << "    -nonorm	 	Read coverage will not be normalized based on sequencing coverage of non-peak loci." << std::endl;
	std::cout << "    -prenorm	 	User must input a pre-normalized data matrix of sequencing coverage as specified in manual." << std::endl;
	std::cout << "    -proc	 	Explicitly state number of processors to use." << std::endl;
	std::cout << "    -lin	 	Perform linear regression." << std::endl;

	std::cout << "\nDetailed usage:\n" << std::endl;
	std::cout << "    Version flag: -v " << std::endl;
	std::cout << "		Display program version and exit program. No argument required.\n" << std::endl; 
	std::cout << "    Help flag: -h" << std::endl;
	std::cout << "		Display this page and exit program. No argument required.\n" << std::endl; 
	std::cout << "    Experiment flag: -bam <directory> " << std::endl;
	std::cout << "		REQUIRED. User specified folder containing sorted bam turnover files including index files." << std::endl;
	std::cout << "		Bam files must be named according with a terminal \"_XXX.bam\", where XXX is an integer time point." << std::endl;
	std::cout << "		If replicates are available, they must contain the same time points. Run replicates separately if" << std::endl;
	std::cout << "		this is not the case\n" << std::endl;
	std::cout << "    Input flag: -i <directory>" << std::endl;	
	std::cout << "		User specified folder containing sorted bam turnover files including index files." << std::endl;
	std::cout << "		Bam files must be named according with a terminal \"_XXX.bam\", where XXX is an integer time point." << std::endl;
	std::cout << "		The time points must match those in the experimental directory. If replicates are available," << std::endl;
 	std::cout << "		ensure the order that input folders are specified are the same order that experimental folders are specified.\n" << std::endl;
	std::cout << "    Bed file flag: -bed <text_file>" << std::endl; 
	std::cout << "		REQUIRED. User specified bed file containing loci of interest. Each line in the file must three tab delimited" << std::endl;
	std::cout << "		columns specifying chromosome, start position and end position. Ensure the bed file is sorted, such as:" << std::endl;
	std::cout << "		sort -k1,1 -k2,2n unsorted.bed > sorted.bed\n" << std::endl;
	std::cout << "    Genome flag: -g <string>" << std::endl;
	std::cout << "		Genome name. Currently supported genomes are human (hg19, hg38), mouse (mm9, mm10), fly (dm3, dm6)," << std::endl;
	std::cout << "		worm (ce11), and yeast (sacCer3).\n" << std::endl;
	std::cout << "    Gene specific 3D profile flag: -3d <text_file>" << std::endl;
	std::cout << "		User specified gene file containing RefSeq gene names. Each gene must be on a separate line\n" << std::endl;
	std::cout << "    Depth flag: -s <double>"  << std::endl;    
	std::cout << "		Saturation threshold consideration for data modelling (allowable range from 0.5-0.95). Default = 0.85." << std::endl;
	std::cout << "		This parameter is important in the clasification of loci type (hills, rises, valleys, and falls). Consequently," << std::endl;
	std::cout << "		this flag is only applicable when the - model flag is set.\n" << std::endl;
	std::cout << "    Leading/trailing flag: -t <int>"  << std::endl;          
	std::cout << "		Time point threshold consideration for data modelling (allowable range from 0-2). Default = 1." << std::endl;
	std::cout << "		This parameter is important in the clasification of loci type (hills, rises, valleys, and falls). Consequently," << std::endl;
	std::cout << "		this flag is only applicable when the - model flag is set.\n" << std::endl;
	std::cout << "    Name flag: -name <string>" << std::endl;
	std::cout << "		User specified name for output files. if not specified, output files will be named: turnover.exp\n" << std::endl;
	std::cout << "    Model flag: -model "  << std::endl;
	std::cout << "		Data modelled based on prediction. Default is no data modelling. If set TDCA predicts loci to be hills, rises" << std::endl;
	std::cout << "		valleys, or falls. If not set, loci are modeled as a rise or fall only (single 5 parameter sigmoid).\n" << std::endl;
	std::cout << "    Depth matrix flag: -dm <text_file> "  << std::endl;
	std::cout << "		Data matrix used to normalize user defined input files. This file must contain integers equal to the number" << std::endl;
	std::cout << "		of total bam files. The reported integers are then applied as total coverage values for each bam file in order" << std::endl;
	std::cout << "		of replicate and then in chronological order.\n" << std::endl;
	std::cout << "    No normalization flag: -nonorm "  << std::endl;
	std::cout << "		Coverage values at each locus will be modelled as is. Normalization based on the sequencing coverage of non-peak " << std::endl;
	std::cout << "		loci will not be performed.\n" << std::endl;
	std::cout << "    Pre-normalized data flag: -prenorm <text_file> "  << std::endl;
	std::cout << "		A text file of loci with additional columns specifying time in tab delimited format. Sequencing coverage data will " << std::endl;
	std::cout << "		be modelled as is.\n" << std::endl;
	std::cout << "    Linear regression flag: -lin"  << std::endl;
	std::cout << "		TDCA will perform linear regression at each locus using all time points and will output relevant parameters." << std::endl;
	return 0;
}
