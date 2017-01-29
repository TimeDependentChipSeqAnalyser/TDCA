/**
 **  This file handles all flags and records necessary values for downstream calculations
 **/
#include <fstream>     		// write to file
#include <iostream>    		// input/output
#include <sstream>
#include <string>
#include <dirent.h>    		// directory access
#include <vector>
#include "LocusTurnover.h"
#include <algorithm>
#include <math.h>

// Forward declarations
int welcome();
int help();
int fileCount(std::string s_dir);
int bamFileCount(std::string s_dir, int i_fileNumber);
int timeLogger(int *arr, std::string s_dir, int i_fileNumber);
double samDregion(std::string s_Bamfile, std::string s_PeakCoordinates);
std::string bedToCoordinates(std::string bedLine);
std::string bamGetter(std::string s_dir, int i_turnoverTime);
double samDtotal(std::string s_Bamfile);
std::string getExecPath();
int peakCounter(std::string s_bed);
int writeReadsToVectorParallel(std::string bamFolderArr[], std::string inputBamFolderArr[], bool b_bamInputFlagCall, int i_bamRep, int i_iterator, std::string s_name, int i_bamFiles, int turnoverTimes[], std::string s_bed, std::vector<std::vector<std::string> > &normReadsVector, bool b_genelistFlagCall, std::string s_genelist, std::vector<double> &normGeneVector, std::string s_genelistFile, std::vector<std::string> &validGeneVector, int i_peakNumber, bool b_norm, double depthMatrix[], bool b_dmFlagCall);
void normVectorAverage(int i_peakNumber, int i_bamFiles, int i_bamRep, std::vector<std::vector<std::string> > &normReadsVector, std::vector<std::vector<std::string> > &normAveVector);
void normVectorError(int i_peakNumber, int i_bamFiles, int i_bamRep, std::vector<std::vector<std::string> > &normReadsVector, std::vector<std::vector<std::string> > &normAveVector, std::vector<std::vector<std::string> > &normErrorVector);
bool genomeGetter(std::string s_dir, std::string s_genome, int i_fileNumber);
int extFileCount(std::string s_dir, int i_fileNumber, std::string s_extension);
std::string bedGetter(std::string s_dir, int i_iterator);
void bedToolsIntersect(std::string s_genomePath, std::string s_bed, std::vector<std::string> &featureVector, int i_loop);
void drcScriptGenerator(int i_loopNumber, int i_bamFiles, int i_truncDiff, int turnoverTimes[], std::vector<double> &turnoverArr, bool b_real, bool b_turnoverDirection, double d_overallMax);
void inflectionWriter(std::vector<std::string> &drcInflectionVector, std::vector<std::string> &drcResidualsVector,  std::vector<std::string> &drcUpperAsymVector, std::vector<std::string> &linResVector);
void rBar(std::string s_rPltsName, int i_bamFiles, int turnoverTimes[], int maxIndexArray[], int minIndexArray[], int i_peakNumber);
void rEnd(std::string s_rPltsName,std::string s_rPltsPDF, bool b_genomeFlagCall, int i_bamRep, bool b_model, int i_increase, int i_decrease);
void rBegin(std::string s_rPltsName, std::string s_rPltsPDF, std::string s_args);
void rDensity(std::string s_rPltsName, int i_peakNumber, std::vector<std::string> drcInflectionVector, int i_bamFiles, int turnoverTimes[], std::string s_name);
void drcNormScriptGenerator(int i_loopNumber, int i_bamFiles, int i_truncDiff, int turnoverTimes[], std::vector<double> &turnoverArr, bool b_real, int i_minIndex, bool b_turnoverDirection, double d_overallMax);
void drcNormScriptGeneratorString(int i_loopNumber, int i_bamFiles, int i_truncDiff, int turnoverTimes[], std::vector<std::vector<std::string> > &normReadsVector, bool b_real, int i_minIndex, int i_bamRep, int i_peakNumber);
void rBoxplot(int turnoverTimes[], std::string s_rPltsName, int i_peakNumber, int i_genomeFileCount, int i_bamFiles, std::vector<std::vector<std::string> > &dataArray, std::string s_direction);
void rIdeohg19(std::string s_rPltsName, int i_bamFiles, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, std::string s_direction);
void rIdeohg38(std::string s_rPltsName, int i_bamFiles, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, std::string s_direction);
void rIdeomm10(std::string s_rPltsName, int i_bamFiles, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, std::string s_direction);
void rIdeomm9(std::string s_rPltsName, int i_bamFiles, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, std::string s_direction);
void rIdeodm3(std::string s_rPltsName, int i_bamFiles, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, std::string s_direction);
void rIdeodm6(std::string s_rPltsName, int i_bamFiles, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, std::string s_direction);
void rIdeoce11(std::string s_rPltsName, int i_bamFiles, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, std::string s_direction);
void normGeneVectorAverage(int i_bamFiles, int i_bamRep, int i_validgenes, std::vector<double> &normGeneVector, std::vector<double> &normGeneAveVector);
void r3Dadjusted(std::string s_genelistFile, std::string s_rGenePltsName, int i_bamFiles, int turnoverTimes[], int i_validgenes, std::vector<double> &normGeneVector, std::vector<std::string> &validGeneVector);
void geneVectorShrinkage(int i_validgenes, int i_bamFiles, std::vector<double> &normGeneVector, std::vector<double> &shrunkGeneVector);
void rBoxRes(std::string s_rPltsName, int i_peakNumber, int i_bamFiles, std::vector<std::string> &truncLinResVector, std::vector<std::string> &linResVector, std::vector<std::vector<std::string> > &dataArray);
void rPie(std::string s_rPltsName, int i_peakNumber, int i_undef, int i_elim, int i_valley, int i_hill, int i_fall, int i_rise, bool b_model);
void rStdev(std::string s_rPltsName, int i_bamFiles, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, int turnoverTimes[], int i_undefRise, int i_undefFall, int i_valley, int i_hill, int i_fall, int i_rise);
void rScat(std::string s_rPltsName, int i_bamFiles, int turnoverTimes[], int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, int i_undefRise, int i_undefFall, int i_valley, int i_hill, int i_fall, int i_rise, bool b_model, std::string s_direction);
void inflectionWriterReverse(std::vector<std::string> &drcInflectionVector, std::vector<std::string> &drcResidualsVector, std::vector<std::string> &linResVector);
std::string exec(const char* cmd);
void drcVersitile(int i_drcPar, int i_maxThreads, int i_loopNumber, int i_bamFiles, int turnoverTimes[], std::vector<double> &turnoverArr, std::vector<bool> &b_typeArr, int i_minIndex, int i_maxIndex, double d_satthresh, std::vector<bool> &b_forwardArr, std::vector<bool> &b_reverseArr, std::vector<bool> &b_peak1Arr, std::vector<bool> &b_peak2Arr, std::vector<bool> &b_dip1Arr, std::vector<bool> &b_dip2Arr, bool b_model);
void inflectionWriterVersitile(int i_maxThreads, int i_peakNumber, std::vector<std::string> &s_drcFI_vec, std::vector<std::string> &s_drcFH_vec, std::vector<std::string> &s_drcFU_vec, std::vector<std::string> &s_drcFL_vec, std::vector<std::string> &s_drcFE_vec, std::vector<std::string> &s_drcRI_vec, std::vector<std::string> &s_drcRH_vec, std::vector<std::string> &s_drcRU_vec, std::vector<std::string> &s_drcRL_vec, std::vector<std::string> &s_drcRE_vec, std::vector<std::string> &s_type_vec);
void rUpperA(std::string s_rPltsName, int i_peakNumber, int i_bamFiles, std::vector<std::vector<std::string> > &dataArray, int i_valley, int i_hill);
void rProfiles(std::string s_rPltsName, int i_peakNumber, int i_bamFiles, std::vector<std::vector<std::string> > &dataArray, bool b_model, double d_satthresh, int turnoverTimes[], int i_valley, int i_hill, int i_fall, int i_rise, int i_undefRise, int i_undefFall);
void rIdeosacCer3(std::string s_rPltsName, int i_bamFiles, int i_peakNumber, std::vector<std::vector<std::string> > &dataArray, std::string s_direction);
std::string drcR(int i_loopNumber, std::string s_peakExtention);
void pkgCheck();
void rHeatmap(std::string s_rPltsName, int i_peakNumber, int i_bamFiles, std::vector<std::vector<std::string> > &dataArray, int turnoverTimes[]);

bool exists(std::string filename)
{

	std::string s_command = "if [ -e ./" + filename + " ]; then echo \"t\"; fi;";
	std::string s_out = exec(s_command.c_str());
	s_out.erase(std::remove(s_out.begin(), s_out.end(), '\n'), s_out.end());
	if (s_out=="t") {
		return true;
	} else {return false;}
}

int main(int argc, char* argv[])
{

	/** -b flag variables **/
	std::string s_bamFlag = "-bam";	//name flag
	bool b_bamFlagCall(false);		//true if bam flag called
	std::string s_bamFolder; 		//folder name containing bam files of different 
							//turnover experiments 
	int i_totalBamFiles; 			//total files in the bamFolder
	int i_bamRep = 0;				//Number of replicates determined from number 
							//of bamFolders
	int i_bamFiles;				//number of bam files

	/** -i flag variables **/
	std::string s_bamInputFlag = "-i";	//name flag
	bool b_bamInputFlagCall(false);	//true if bamInput flag called
	std::string s_inputBamFolder; 	//folder name containing input bam files
							//for different turnover experiments.  
	int i_totalInputBamFiles;		//total files in the inputBamFolders
	int i_bamInputRep = 0; 			//Number of replicates determined from number 
							//of inputBamFolders. bamFolders == inputBamFolders.
	int i_bamInputFiles; 			//number of input bam files

	/** -g flag variables **/
	std::string s_genomeFlag = "-g";	//name flag
	bool b_genomeFlagCall(false);		//true if genome flag called
	std::string s_genome; 			//what genome was called

	/** -3d flag variables **/
	std::string s_genelistFlag = "-3d";	//name flag
	bool b_genelistFlagCall(false);	//true if gene list flag called
	std::string s_genelist; 		//file name of gene list

	/** -bed flag variables **/
	std::string s_bedFlag = "-bed";	//name flag
	bool b_bedFlagCall(false);		//true if bed flag called (if false, end program)
	std::string s_bed; 			//name of bed file

	/** -name, --name flag variables **/
	std::string s_nameFlag = "-name";	//name flag
	std::string s_name; 			//output name. If std::string = NULL, call output 
							//files some default name

	/** -v, --version flag variables **/
	bool b_versionFlagCall(false);	//if true, only print version and do not run 
							//program regardless of other flags

	/** -h, --help flag variables **/
	bool b_helpFlagCall(false);		//if true, only print help page and do not run 
							//program regardless of other flags

	/** -s, --sat flag saturation threshold **/
	double d_satthresh = 0.85;		// default is 0.85
	std::string s_satFlag = "-s";
	bool b_satFlagCall(false);		// false if user has not specified desired value

	/** -model flag  **/
	bool b_model(false);			//if true, pedict how data should be modelled
	std::string s_modelFlag = "-model";

	/** -norm flag  **/
	bool b_norm(true);			//if false, do not normalize based on total depth
	std::string s_normFlag = "-nonorm";

	/** -trail flag  **/
	std::string s_trailFlag = "-t"; //if not specified, i_trailCutoff = 1
	int i_trailCutoff = 1;
	bool b_trailFlagCall = false;

	/** -dm flag  **/
	std::string s_dmFlag = "-dm";
	std::string s_dmFile;
	bool b_dmFlagCall = false;

	/** -p flag  **/
	std::string s_pFlag = "-p";
	int i_threads = 1;
	bool b_pFlagCall = false;

	std::string s_execPath = getExecPath();
	std::cout << "Program path: " << s_execPath << std::endl;

	// display welcome message when there is no argument given
	if (argc < 2) {
		welcome();
		std::exit(0); }

	// count bam folder and input bam folder reps
	for (int i = 1; i < argc; i++){
		std::string arg = argv[i];
		if (arg == "-bam"){
			b_bamFlagCall = true;
			i_bamRep++; }
		if (arg == "-i"){
			b_bamInputFlagCall = true;
			i_bamInputRep++; }
	}

	// make an array of all bam folder names
	std::string bamFolderArr[i_bamRep];
	int i_bamFolderPos = 0; //index place holder for bamFolderArr

	std::string inputBamFolderArr[i_bamInputRep];
	int i_inputBamFolderPos = 0; //index place holder for inputBamFolderArr 

	// Go through arguments again and check if input is specified after appropriate flags
	for (int i = 1; i < argc; i++){
		// print out each argument
		std::string arg = argv[i];

		if ((arg == "--help") || (arg == "-h") || (arg == "-help")){
			b_helpFlagCall = true;
		} else if ((arg == "--version") || (arg == "-v") || (arg == "-version")) {
			b_versionFlagCall = true;	
		} else if (arg.find(s_bamFlag) != std::string::npos) {
			i++;
			try {	std::string s = argv[i]; }
			catch(...){
				std::cerr << "-bam flag must preceed a folder name" << std::endl;
				std::exit(0); }
			bamFolderArr[i_bamFolderPos] = argv[i];
			i_bamFolderPos++;
		} else if (arg.find(s_bamInputFlag) != std::string::npos) {
			i++;
			try {	
				std::string s = argv[i]; }
			catch(...){
				std::cerr << "-i flag must preceed a folder name" << std::endl;
				std::exit(0); }
			inputBamFolderArr[i_inputBamFolderPos] = argv[i];
			i_inputBamFolderPos++;
		} else if (arg.find(s_bedFlag) != std::string::npos) {
			i++;
			b_bedFlagCall = true;
			try {	
				std::string s = argv[i];
				bool b_ex = exists(s);
				if (b_ex == 0) {
					std::cout << s << " is not a file in directory." << std::endl;
					std::exit(0); }
			}
			catch(...){
				std::cerr << "-bed flag must preceed a file name" << std::endl;
				std::exit(0); }
			s_bed = argv[i];
		} else if (arg.find(s_genomeFlag) != std::string::npos) {
			i++;
			b_genomeFlagCall = true;
			try {	std::string s = argv[i]; }
			catch(...){
				std::cerr << "-g flag must preceed a genome name" << std::endl;
				std::exit(0); }
			s_genome = argv[i];

			// check if genome input is valid
			if(genomeGetter(s_execPath + "lib/GenomeFeatures/", s_genome, fileCount(s_execPath + "lib/GenomeFeatures/"))) {
				std::cout << "Valid genome feature folder detected." << std::endl;
			} else {
				std::cerr << s_genome << " is not a valid genome." << std::endl;
				std::exit(0); }
		} else if (arg.find(s_nameFlag) != std::string::npos) {
			i++;
			try {	std::string s = argv[i]; }
			catch(...){
				std::cerr << "-n flag must preceed an output file name" << std::endl;
				std::exit(0); }
			s_name = argv[i];
		} else if (arg.find(s_genelistFlag) != std::string::npos){
			i++;
			b_genelistFlagCall = true;
			try {	
				std::string s = argv[i];
				bool b_ex = exists(s);
				if (b_ex == 0) {
					std::cout << s << " is not a file in directory." << std::endl;
					std::exit(0); }
			}
			catch(...){
				std::cerr << "-3d flag must preceed a file name" << std::endl;
				std::exit(0); }
			s_genelist = argv[i];
		} else if (arg.find(s_trailFlag) != std::string::npos) {
			i++;
			try {	
				std::string s = argv[i];
				i_trailCutoff = std::stoi(argv[i]);
				b_trailFlagCall = true;
				if (i_trailCutoff > 2 || d_satthresh < 0) {
					std::cerr << "-t flag must be between 0 and 2." << std::endl;
					std::exit(0); }
			}
			catch(...){
				std::cerr << "-t flag must be between 0 and 2." << std::endl;
				std::exit(0); }
		} else if (arg.find(s_satFlag) != std::string::npos) {
			i++;
			try {	
				std::string s = argv[i];
				d_satthresh = std::stod(argv[i]);
				b_satFlagCall = true;
				if (d_satthresh > 1 || d_satthresh < 0.5) {
					std::cerr << "-s flag must be between 0.95 and 0.5." << std::endl;
					std::exit(0); }
			}
			catch(...){
				std::cerr << "-s flag must be between 0.95 and 0.5." << std::endl;
				std::exit(0); }
		} else if (arg == s_modelFlag) {
			b_model = true;
		} else if (arg.find(s_dmFlag) != std::string::npos){
			i++;
			b_dmFlagCall = true;
			try {	
				std::string s = argv[i];
				bool b_ex = exists(s);
				if (b_ex == 0) {
					std::cout << s << " is not a file in directory." << std::endl;
					std::exit(0); }
			}
			catch(...){
				std::cerr << "-3d flag must preceed a file name" << std::endl;
				std::exit(0); }
			s_dmFile = argv[i];
		} else if (arg == s_normFlag) {
			b_norm = false;
		} else if (arg.find(s_pFlag) != std::string::npos) {
			i++;
			try {	
				std::string s = argv[i];
				i_threads = std::stoi(argv[i]);
				b_pFlagCall = true;
				if (i_threads < 1 || i_threads > 100) {
					std::cerr << "-p flag must be an integer between 1 and 100." << std::endl;
					std::exit(0); }
			}
			catch(...){
				std::cerr << "-p flag must be an integer between 1 and 100." << std::endl;
				std::exit(0); }
		} else {
			std::cout << "Invalid flag called. Program terminated." <<  std::endl;
			std::exit(0); }
	} //for loop

	// check if user specified information is logical
	if (b_helpFlagCall) {
		help();
		std::exit(0); }
	if (b_versionFlagCall) {
		std::cout << "tdca 1.1.0.01-28-2017" <<  std::endl;
		std::exit(0); }
	if (!b_bedFlagCall) {
		std::cout << "Bed file required. Program terminated." <<  std::endl;
		std::exit(0); }
	if (!b_bamFlagCall) {
		std::cout << "Bam folder required. Program terminated." <<  std::endl;
		std::exit(0); }
	if (i_bamInputRep > 0 && i_bamRep != i_bamInputRep) {
		std::cout << "Bam replicates do not equal bam input replicates. Program terminated." << std::endl;
		std::exit(0); }
	if (b_genelistFlagCall == true && b_genomeFlagCall == false) {
		std::cout << "Gene list flag (-3d) must be combined with genome flag (-g). Program terminated." << std::endl;
		std::exit(0); }
	if (b_trailFlagCall == true && b_model == false) {
		std::cout << "Trailing does not apply to data that is not modelled (cannot run -t without -model). Program terminated." << std::endl;
		std::exit(0); }
	if (b_satFlagCall == true && b_model == false) {
		std::cout << "Saturation does not apply to data that is not modelled (cannot run -s without -model). Program terminated." << std::endl;
		std::exit(0); }



	// bam folder checks to ensure number of bam files are equivalent
	for (int i = 0; i < i_bamRep; i++) {

		std::string s_relativePath ("./");
		std::string s_bamDirectory1 = s_relativePath + bamFolderArr[i];

		//number of bam files in directory
		int i_bamFiles1 = bamFileCount(s_bamDirectory1, fileCount(s_bamDirectory1));

		if (i_bamFiles1 < 4) {
			std::cout << "Turnover values cannot be calculated with less" <<
				"then four time points." << std::endl;
			std::exit(0); }

		// check if bam folder replicates have same number of bam files
		if (i < i_bamRep-1) {
			std::string s_bamDirectory2 = s_relativePath + bamFolderArr[i+1];
			int i_bamFiles2 = bamFileCount(s_bamDirectory2, fileCount(s_bamDirectory2));

			if (i_bamFiles1 != i_bamFiles2) {
				std::cout << "Number of bam files not equivalent in each bam" <<
					" folder replicate." << std::endl;
				std::exit(0); }
		}

		// check if bam folder has same number of bam files in input if available
		if (i_bamInputRep > 0) {
			std::string s_inputBamDirectory = s_relativePath + inputBamFolderArr[i];
			int i_inputBamFiles = bamFileCount(s_inputBamDirectory, fileCount(s_inputBamDirectory));

			if (i_bamFiles1 != i_inputBamFiles) {
				std::cout << "Number of bam files not equivalent in each bam" <<
					" folder and bam input folder." << std::endl;
				std::exit(0); }
		}
	}

	// check if time extensions are equal in each bam folder
	std::string s_relativePath ("./");
	std::string s_bamDirectory = s_relativePath + bamFolderArr[0];
	i_bamFiles = bamFileCount(s_bamDirectory, fileCount(s_bamDirectory));

	int turnoverTimes[i_bamFiles]; // array containing ChIP-seq turnover times

	timeLogger(turnoverTimes, s_bamDirectory, fileCount(s_bamDirectory));

	// sort the turnoverTimes array in asending order
	int temp = 0;
	for(int i = 0; i < i_bamFiles; i++) {
		for(int j = 0; j < i_bamFiles; j++) {
			if(turnoverTimes[i] < turnoverTimes[j]) {
				temp = turnoverTimes[i];
				turnoverTimes[i] = turnoverTimes[j];
				turnoverTimes[j] = temp; }
		}
	}

	for (int i = 0; i < i_bamRep; i++) {

		// check if time extensions are equal in each bam folder
		std::string s_relativePath ("./");
		std::string s_bamDirectory = s_relativePath + bamFolderArr[i];
		i_bamFiles = bamFileCount(s_bamDirectory, fileCount(s_bamDirectory));

		int otherTurnoverTimes[i_bamFiles]; // array containing ChIP-seq turnover 
						    // times for all other bam folders
		timeLogger(otherTurnoverTimes, s_bamDirectory, fileCount(s_bamDirectory));

		// sort the otherTurnoverTimes array in asending order
		int temp = 0;
		for(int i = 0; i< i_bamFiles; i++) {
			for(int j = 0; j < i_bamFiles; j++) {
				if(otherTurnoverTimes[i] < otherTurnoverTimes[j]) {
					temp = otherTurnoverTimes[i];
					otherTurnoverTimes[i] = otherTurnoverTimes[j];
					otherTurnoverTimes[j] = temp; }
			}
		}

		for (int i = 0; i < i_bamFiles; i++) {
			if (turnoverTimes[i] != otherTurnoverTimes[i]) {
				std::cout << "Bam turnover times not equivalent in each bam" <<
						" folder replicate." << std::endl;
				std::exit(0); }
		}

		// check if bam folder has same turnover extensions as input if available
		if (i_bamInputRep > 0) {
			std::string s_bamDirectory = s_relativePath + inputBamFolderArr[i];
			i_bamFiles = bamFileCount(s_bamDirectory, fileCount(s_bamDirectory));

			int otherTurnoverTimes[i_bamFiles]; // array containing ChIP-seq turnover 
							    // times for all other bam foders
			timeLogger(otherTurnoverTimes, s_bamDirectory, fileCount(s_bamDirectory));

			// sort the otherTurnoverTimes array in asending order
			int temp = 0;
			for(int i = 0; i< i_bamFiles; i++) {
				for(int j = 0; j < i_bamFiles; j++) {
					if(otherTurnoverTimes[i] < otherTurnoverTimes[j]) {
						temp = otherTurnoverTimes[i];
						otherTurnoverTimes[i] = otherTurnoverTimes[j];
						otherTurnoverTimes[j] = temp; }
				}
			}

			for (int i = 0; i < i_bamFiles; i++) {
				if (turnoverTimes[i] != otherTurnoverTimes[i]) {
					std::cout << "Bam turnover times not equivalent in bam experiment" <<
							" folder and bam input folder." << std::endl;
					std::exit(0); }
			}
		}
	} 

	// check if drc is installed
	pkgCheck();

	std::string s_args = "Bam rep = "; // write to R script at begining
	// output all user specified information
	std::cout << "TDCA will run based on the following parameters:" << std::endl;

	if (b_bamFlagCall) {
		std::cout << "Bam replicates = " << i_bamRep << std::endl; 
		for (int i = 0; i < i_bamRep; i++) {
			std::cout << "Bam rep " << i+1 << " = " << bamFolderArr[i] << std::endl;
			s_args += bamFolderArr[i] + ", "; }
	}
	if (b_bamInputFlagCall) {
		std::cout << "Input bam replicates = " << i_bamInputRep << std::endl; 
		s_args += "Bam input rep = ";
		for (int i = 0; i < i_bamInputRep; i++) {
			std::cout << "Input bam rep " << i+1 << " = " << inputBamFolderArr[i] << std::endl;
			s_args += inputBamFolderArr[i] + ", "; }
	}
	if (b_bedFlagCall) {
		std::cout << "Bed file name = " << s_bed << std::endl;
		s_args += "Bed file name = " + s_bed + ", "; }
	if (b_genomeFlagCall) {
		std::cout << "Genome name = " << s_genome << std::endl;
		s_args += "Genome name = " + s_genome + ", "; }
	if (!s_name.empty()) {
		std::cout << "Output files will be named " << s_name << std::endl;
		s_args += "Output files will be named " + s_name + ", ";
	} else {
		s_name = "turnover.exp";
		std::cout << "Output files will be named " << s_name << " by default" << std::endl; }
	if (b_model) {
		std::cout << "Data will be modelled based on predicted behaviour." << std::endl; 
		std::cout << "-t parameter set as " << i_trailCutoff << std::endl; }
	if ( (b_model == true) && (b_satFlagCall == false) ) {
		std::cout << "-s parameter set as "<< d_satthresh << std::endl; }
	if (b_satFlagCall) {
		std::cout << "-s parameter set as "<< d_satthresh << std::endl; }
	if (b_dmFlagCall) {
		std::cout << "-dm flag called. TDCA will use depth matrix file for total depth values." << std::endl; }
	if (b_pFlagCall) {
		std::cout << "-p flag called. TDCA will use " << i_threads << " for parallelization." << std::endl; }

	std::cout << "The following turnover times were detected: " << std::endl;
	s_args += "The following turnover times were detected: ";
	for (int i = 0; i < i_bamFiles; i++) { // for each bam file
		std::cout << "time " << i+1 << " = " << turnoverTimes[i] << "." << std::endl;
		s_args += std::to_string(turnoverTimes[i]) + " "; }

	// print number of peaks in bed file
	int i_peakNumber = peakCounter(s_bed);
	std::cout << "The bed file " << s_bed << " contains " << i_peakNumber << " peaks." << std::endl;

	// make an array of files in genome folder for gene feature overlap
	std::string s_genomeDir = s_execPath + "lib/GenomeFeatures/" + s_genome + "/";
	// XXX THUMBNAIL FILES ARE COUNTED XXX //


	int i_genomeFileCount = extFileCount(s_genomeDir, fileCount(s_genomeDir), ".bed");
	std::string genomeFilesArr[i_genomeFileCount];
	std::cout << i_genomeFileCount << " gene feature files will be used:" <<std::endl;	
	for (int i = 0; i < i_genomeFileCount; i++) {
		genomeFilesArr[i] = bedGetter(s_genomeDir, i);
		std::cout << genomeFilesArr[i] << std::endl; }


 	// ======================== XXX CALCULATION COMPONENT OF PROGRAM XXX =======================//

	// type vectors
	std::vector <bool> b_forwardArr = {0,0,0,1,0,0}; 
	std::vector <bool> b_reverseArr = {0,1,0,0,0,0}; 
	std::vector <bool> b_peak1Arr   = {0,1,0,0,1,0};
	std::vector <bool> b_peak2Arr   = {0,0,1,1,0,0};  
	std::vector <bool> b_dip1Arr    = {1,0,0,1,0,0}; 
	std::vector <bool> b_dip2Arr    = {0,1,0,0,0,1}; 

	//================================== XXX VARIABLE RECAP XXX ================================//
	// Number of bam file turnover times = i_bamFiles
	// Turnover times stored in ordered array turnoverTimes[]. Use this as REGEX to samD bam files
	// Number of replicates = i_bamRep
	// Array storing folder name of replicates = bamFolderArr[]
	// Array storing folder name of replicates = inputBamFolderArr[]
	//==========================================================================================//

	// vector holding ALL loci information - print contents to output file.
	// rows = i_peakNumber + 1; columns	= 2*i_bamFiles + i_genomeFileCount + 11
	int i_dataArrCols = (2*i_bamFiles + i_genomeFileCount + 12);
	std::vector< std::vector<std::string> > dataArray ((i_peakNumber + 1), std::vector<std::string>(i_dataArrCols, "-" ) );

	// make the header for dataArray 	
	dataArray[0][0] = "chrom";							// chromosome
	dataArray[0][1] = "start";							// start
	dataArray[0][2] = "end";							// end
	for (int i = 3; i < (i_bamFiles + 3); i++) 
		dataArray[0][i] = ("time = " + std::to_string(turnoverTimes[i-3]));

	for (int i = (3 + i_bamFiles); i < (2*i_bamFiles + 3); i++) 
		dataArray[0][i] = ("time = " + std::to_string(turnoverTimes[i-(3 + i_bamFiles)]) + " stdev");

	dataArray[0][2*i_bamFiles + 3]  = "s_type_vec[i]"; 
	dataArray[0][2*i_bamFiles + 4]  = "s_drcFI_vec[i]"; 
	dataArray[0][2*i_bamFiles + 5]  = "s_drcFH_vec[i]";
	dataArray[0][2*i_bamFiles + 6]  = "s_drcFU_vec[i]"; 
	dataArray[0][2*i_bamFiles + 7]  = "s_drcFL_vec[i]";	
	dataArray[0][2*i_bamFiles + 8]  = "s_drcRI_vec[i]"; 
	dataArray[0][2*i_bamFiles + 9]  = "s_drcRH_vec[i]"; 
	dataArray[0][2*i_bamFiles + 10] = "s_drcRU_vec[i]"; 
	dataArray[0][2*i_bamFiles + 11] = "s_drcRL_vec[i]"; 
	
	for ( int i = (i_dataArrCols - i_genomeFileCount); i < (i_dataArrCols); i++)
		dataArray[0][i] = genomeFilesArr[i-(i_dataArrCols - i_genomeFileCount)]; // feature hits

	// Bedtools intesect on genome features commands - make additional bed files of genome features
	if (b_genomeFlagCall) {
		std::cout << "Intersecting peaks with genome features." << std::endl;
		#pragma omp parallel for 
		for (int i = 0; i < i_genomeFileCount; i++) {
			std::vector<std::string> featureVector(i_peakNumber, "0");
			//std::cout << "s_genomeDir + genomeFilesArr[i] = " << s_genomeDir + genomeFilesArr[i] << std::endl;
			bedToolsIntersect(s_genomeDir + genomeFilesArr[i], s_bed, featureVector, i);

			for (int j = 0; j < i_peakNumber; j++)
				dataArray[(j+1)][(i_dataArrCols - i_genomeFileCount + i)] = featureVector[j]; 
		}
	}

	// write normalized reads to vector
	std::vector<std::vector<std::string> > normReadsVector;

	// write normalized gene reads to vector
	std::vector<double> normGeneVector;

	// write the header to normReadsVector
	std::vector<std::string> headerVector;
	headerVector.push_back("Region");

	for (int j = 0; j < i_bamFiles; j++) // for each bam file
		headerVector.push_back(std::to_string(turnoverTimes[j]));
	normReadsVector.push_back(headerVector);

	// establish directory of proper genome refseq gene list
	std::string s_genelistFile = s_execPath + "lib/GeneLists/" + s_genome + "/" + s_genome + ".genes";
	// vector containing VALID user input refseq gene names
	std::vector<std::string> validGeneVector;
	// number of user input VALID refseq gene names
	int i_validgenes;

	// dm matrix
	int i_inputMultiplier;
	if (b_bamInputFlagCall) { i_inputMultiplier = 2 ; }
	else  { i_inputMultiplier = 1 ; }

	double depthMatrix[i_bamFiles*i_inputMultiplier*i_bamRep];
	if (b_dmFlagCall) {
		std::ifstream dmFile;
		dmFile.open(s_dmFile);
		std::string s_line;
		for (int i = 0; i < i_bamFiles*i_inputMultiplier*i_bamRep; i++) {
			std::getline(dmFile,s_line);
			try {
				depthMatrix[i] = stod(s_line);
			}
			catch(...){
				std::cerr << "Values in depth matrix file must be integers. Program terminated" << std::endl;
				std::exit(0);
			}
			if (!(depthMatrix[i] > 0)) {
				std::cerr << "Values in depth matrix file must be greater than 0. Program terminated" << std::endl;
				std::exit(0);
			}
		}
		dmFile.close();
	}

	for (int i = 0; i < i_bamRep; i++) // for each folder
		i_validgenes = writeReadsToVectorParallel(bamFolderArr, inputBamFolderArr, b_bamInputFlagCall, i_bamRep, i, s_name, i_bamFiles, turnoverTimes, s_bed, normReadsVector, b_genelistFlagCall, s_genelist, normGeneVector, s_genelistFile, validGeneVector, i_peakNumber, b_norm, depthMatrix, b_dmFlagCall);

	// Compress normalized reads in normvector if there is mltiple replicates
	std::vector<std::vector<std::string> > normAveVector;

	// Create a vector to hold variability of compressed normalized reads when multiple replicates are available
	std::vector<std::vector<std::string> > normErrorVector;

	// compress replicates into one vector and give an additional vector for replicate variance
	if (i_bamRep > 1) {
		//std::cout << "Compressing normReadsVector to normAveVector." << std::endl;
		normVectorAverage(i_peakNumber, i_bamFiles, i_bamRep, normReadsVector, normAveVector);
		normVectorError(i_peakNumber, i_bamFiles, i_bamRep, normReadsVector, normAveVector, normErrorVector);
	}	

	// Compress normalized reads in normGeneVector if there is mltiple replicates
	std::vector<double> normGeneAveVector;
	if (b_genelistFlagCall == true && i_bamRep > 1)
		normGeneVectorAverage(i_bamFiles, i_bamRep, i_validgenes, normGeneVector, normGeneAveVector);

	// write normReadsVector OR normAveVector && normErrorVector to dataArray
	if (i_bamRep > 1) { // multiple replicates
		for (int i = 1; i < i_peakNumber+1; i++) {
			for (int j = 0; j < i_bamFiles; j++) {
				dataArray[i][j + 3] = normAveVector[i][j+1]; 
				dataArray[i][j + i_bamFiles + 3] = normErrorVector[i][j+1]; }
		}
	} else { // one replicate
		for (int i = 1; i < i_peakNumber+1; i++) {
			for (int j = 0; j < i_bamFiles; j++)
				dataArray[i][j + 3] = normReadsVector[i][j+1];
		}
	}


/**===========================================================================================

ABOVE HERE I BELEIVE EVERYTHING HAS TO DO WITH FLAG HANDLING, DEPTH EXTRACTION< NORMALIZATION AND REPLICATE HANDLING
*EXCEPT* THAT THE NAMES OF THE FIRST ROW IN dataArray MUST BE EDITED APPROPRIATELY

BELOW HERE CHANGES SHOULD BE MADE FOR NEW TDCA 

**///==========================================================================================


	// R script file
	std::string s_rPltsName = s_name + ".tdca.R";
	std::string s_rPltsPDF = s_name + ".tdca.pdf";
	rBegin(s_rPltsName, s_rPltsPDF, s_args); 

	// initialize arrays for histogram of timepoints of min and max 
	int maxIndexArray [i_bamFiles];
	int minIndexArray [i_bamFiles];
	for (int i = 0; i < i_bamFiles; i++) {
		maxIndexArray[i] = 0;
		minIndexArray[i] = 0;	
	}

	// INITIALIZE THESE VECTORS FOR DRC VERSITILE 
	std::vector<std::string> s_type_vec(i_peakNumber);
	std::vector<std::string> s_drcFI_vec(i_peakNumber, "-");
	std::vector<std::string> s_drcFH_vec(i_peakNumber, "-");
	std::vector<std::string> s_drcFU_vec(i_peakNumber, "-");
	std::vector<std::string> s_drcFL_vec(i_peakNumber, "-");
	std::vector<std::string> s_drcFE_vec(i_peakNumber, "-");
	std::vector<std::string> s_drcRI_vec(i_peakNumber, "-");
	std::vector<std::string> s_drcRH_vec(i_peakNumber, "-");
	std::vector<std::string> s_drcRU_vec(i_peakNumber, "-");
	std::vector<std::string> s_drcRL_vec(i_peakNumber, "-");
	std::vector<std::string> s_drcRE_vec(i_peakNumber, "-");


	std::remove("drcVersitile.txt"); 
	std::remove("drcVersitile.R"); 

	int i_maxThreads = 0;
	if (b_pFlagCall) { i_maxThreads = i_threads; }
	else { i_maxThreads = stoi(exec(std::string("nproc").c_str())); }
	int i_drcPar;						 		
/**
	// ===================================================================================================================
	// XXX THIS IS PROBABLY NOT GOOD XXX - maybe make this better in future
	// NEED TO WRITE SEPERATE FUNCTIONS FOR s ESTIMATE. CANOT USE drcVersitile OR inflection WriterVersitile.
	// ESTIMATE SATTHRESH
	// ==================
	if ( (b_satFlagCall == false) && (b_model == true) ) { // skipped if nomodel (no eliminated loci)
		std::cout << "Estimating s value" << std::endl;
		int i_testLoci;
		if (i_peakNumber < 100) { i_testLoci = i_peakNumber; }
		else { i_testLoci = 100; }
		i_drcPar = floor(i_testLoci/i_maxThreads);

		int i_e9 = 0;
		int i_e8 = 0;
		int i_e7 = 0;
		int i_e6 = 0;
		int i_e5 = 0;

		// test the following s values: 0.95, 0.85, 0.75, 0.65, 0.55
		for (int k = 0; k < 5; k++) {
			if (k == 0) { d_satthresh = 0.95; }
			if (k == 1) { d_satthresh = 0.85; }
			if (k == 2) { d_satthresh = 0.75; }
			if (k == 3) { d_satthresh = 0.65; }
			if (k == 4) { d_satthresh = 0.55; }

			int i_loopCounter = 0;
			for (int z = 1; z < i_testLoci+1; z++) {
				int i = (i_peakNumber/i_testLoci)*z;
				std::string objName;
				if (i_bamRep > 1) {
					objName = normAveVector[i][0];
				} else {
					objName = normReadsVector[i][0];
				}
				LocusTurnover ($objName) (objName, i_bamFiles, d_satthresh, i_trailCutoff);
				// write data to objects
				if (i_bamRep > 1) {
					for (int j = 0; j < i_bamFiles; j++) {
						double d_depth = stod(normAveVector[i][j+1]); 
						($objName).turnoverArr[j] = d_depth;
					}
				} else {
					for (int j = 0; j < i_bamFiles; j++) {
						double d_depth = stod(normReadsVector[i][j+1]);
						($objName).turnoverArr[j] = d_depth;
					}
				}
				($objName).analysis();
				int i_maxIndex = ($objName).i_overallMaxIndex;
				int i_minIndex = ($objName).i_overallMinIndex;

				if (($objName).b_typeArr == b_forwardArr ) {
					s_type_vec[i-1] = "rise";
				} else if (($objName).b_typeArr == b_reverseArr) { 
					s_type_vec[i-1] = "fall";
				} else if ( (($objName).b_typeArr == b_peak1Arr) || (($objName).b_typeArr == b_peak2Arr) ) { 
					s_type_vec[i-1] = "hill";
				} else if ( (($objName).b_typeArr == b_dip1Arr) || (($objName).b_typeArr == b_dip2Arr) ) { 
					s_type_vec[i-1] = "valley";
				} else { 
					s_type_vec[i-1] = "undefined"; 
				}

				if (i_bamRep > 500) { // THIS CALCULATES INFLECTION POINT OF REPLICATES SEPERATELY THEN AVERAGES.
					//drcNormScriptGeneratorString(i_loopCounter, i_bamFiles, i_truncDiff, turnoverTimes, normReadsVector, ($objName).b_realPeak, i_minIndex, i_bamRep, i_peakNumber);					
				} else {											
					drcVersitile(i_drcPar, i_maxThreads, i_loopCounter, i_bamFiles, turnoverTimes, ($objName).turnoverArr, ($objName).b_typeArr, i_minIndex, i_maxIndex, d_satthresh, b_forwardArr, b_reverseArr, b_peak1Arr, b_peak2Arr, b_dip1Arr, b_dip2Arr, b_model);
				}	
											
				i_loopCounter++;
			} // end of object loop 
			inflectionWriterVersitile(i_maxThreads, i_testLoci, s_drcFI_vec, s_drcFH_vec, s_drcFU_vec, s_drcFL_vec, s_drcFE_vec, s_drcRI_vec, s_drcRH_vec, s_drcRU_vec, s_drcRL_vec, s_drcRE_vec, s_type_vec);
			// ELIMINATE BAD LOCI
			for (int i = 0; i < i_testLoci; i++) {
				if ( (s_drcFH_vec[i] == "NaN") || (s_drcRH_vec[i] == "NaN") ) {
					if (k == 0) { i_e9++; }
					if (k == 1) { i_e8++; }
					if (k == 2) { i_e7++; }
					if (k == 3) { i_e6++; }
					if (k == 4) { i_e5++; }
				} else if ( ((s_type_vec[i] == "rise") && (stod(s_drcFH_vec[i]) > 0)) || // bad rise
				     ((s_type_vec[i] == "fall") && (stod(s_drcRH_vec[i]) < 0)) || // bad fall
				     ((s_type_vec[i] == "hill") && ((stod(s_drcRH_vec[i]) < 0) || (stod(s_drcFH_vec[i]) > 0))) ||     // bad hill
				     ((s_type_vec[i] == "valley") && ((stod(s_drcRH_vec[i]) < 0) || (stod(s_drcFH_vec[i]) > 0))) ) {  // bad valley
						if (k == 0) { i_e9++; }
						if (k == 1) { i_e8++; }
						if (k == 2) { i_e7++; }
						if (k == 3) { i_e6++; }
						if (k == 4) { i_e5++; }
				} else {}
			}
		}
		std::vector<double> d_eRatioVec;
		d_eRatioVec.push_back(i_e9 * (1/0.95));
		d_eRatioVec.push_back(i_e8 * (1/0.85));
		d_eRatioVec.push_back(i_e7 * (1/0.75));
		d_eRatioVec.push_back(i_e6 * (1/0.65));
		d_eRatioVec.push_back(i_e5 * (1/0.55));

		int i_minIndex = 4;
		//std::cout << "d_eRatioVec[4] = " << d_eRatioVec[4] << std::endl;
		//for (int i = d_eRatioVec.size() - 1; i > -1; i--) {
		//	std::cout << "d_eRatioVec[" << i << "] = " << d_eRatioVec[i] << std::endl;
		//	if (d_eRatioVec[i] < d_eRatioVec[i_minIndex]) { i_minIndex = i; }
		//}
		if (i_minIndex == 0) { d_satthresh = 0.95; }
		if (i_minIndex == 1) { d_satthresh = 0.85; }
		if (i_minIndex == 2) { d_satthresh = 0.75; }
		if (i_minIndex == 3) { d_satthresh = 0.65; }
		if (i_minIndex == 4) { d_satthresh = 0.55; }

		// CLEAN UP VECTORS FOR DRC OBJECT LOOP
		for (int i = 0; i < i_peakNumber; i++) {
			s_type_vec[i] = "-";
			s_drcFI_vec[i] = "-";
			s_drcFH_vec[i] = "-";
			s_drcFU_vec[i] = "-";
			s_drcFL_vec[i] = "-";
			s_drcFE_vec[i] = "-";
			s_drcRI_vec[i] = "-";
			s_drcRH_vec[i] = "-";
			s_drcRU_vec[i] = "-";
			s_drcRL_vec[i] = "-";
			s_drcRE_vec[i] = "-";
		}
	std::cout << "-s parameter estimated as: " << d_satthresh << std::endl;
	} // end of s parameter estimate

	//===================================================================================================================
**/






	int i_loopCounter = 0;
	// for progress
	float f_progress = 0.0;
	int i_barWidth = 20;
	double d_progiterator = 0;

	std::cout << "Writing drc script." << std::endl;

	std::remove("drcVersitile.txt"); 
	std::remove("drcVersitile.R");  

	if (b_pFlagCall) { i_maxThreads = i_threads; }
	else { i_maxThreads = stoi(exec(std::string("nproc").c_str())); }
	i_drcPar = floor(i_peakNumber/i_maxThreads); 
	// OBJECT LOOP
	// ===========
	for (int i = 1; i < i_peakNumber + 1; i++) {
		std::string objName;
		if (i_bamRep > 1) {
			objName = normAveVector[i][0];
		} else {
			objName = normReadsVector[i][0];
		}
		LocusTurnover ($objName) (objName, i_bamFiles, d_satthresh, i_trailCutoff);
		// write data to objects
		if (i_bamRep > 1) {
			for (int j = 0; j < i_bamFiles; j++) {
				double d_depth = stod(normAveVector[i][j+1]); 
				($objName).turnoverArr[j] = d_depth;
			}
		} else {
			for (int j = 0; j < i_bamFiles; j++) {
				double d_depth = stod(normReadsVector[i][j+1]);
				($objName).turnoverArr[j] = d_depth;
			}
		}
		($objName).analysis();
		dataArray[i][0] = ($objName).s_chr;
		dataArray[i][1] = ($objName).s_start;
		dataArray[i][2] = ($objName).s_end;

		if (b_model) {
			if (($objName).b_typeArr == b_forwardArr ) {
				s_type_vec[i-1] = "rise";
			} else if (($objName).b_typeArr == b_reverseArr) { 
				s_type_vec[i-1] = "fall";
			} else if ( (($objName).b_typeArr == b_peak1Arr) || (($objName).b_typeArr == b_peak2Arr) ) { 
				s_type_vec[i-1] = "hill";
			} else if ( (($objName).b_typeArr == b_dip1Arr) || (($objName).b_typeArr == b_dip2Arr) ) { 
				s_type_vec[i-1] = "valley";
			} else { 
				s_type_vec[i-1] = "undefined"; 
			}
		} else {
			s_type_vec[i-1] = "undefined"; // nomodel is undefined so that drcVersitile uses all data. Fix this for output
		}

		// fill in array for histogram of timepoints of min and max 
		int i_maxIndex = ($objName).i_overallMaxIndex;
		for (int m = 0; m < i_bamFiles; m++) {
			if (m == i_maxIndex) {
				int i_temp = maxIndexArray[m]; 
				maxIndexArray[m] = (i_temp + 1);
			}
		} 
		int i_minIndex = ($objName).i_overallMinIndex;
		for (int n = 0; n < i_bamFiles; n++) {
			if (n == i_minIndex) {
				int i_temp = minIndexArray[n];
				minIndexArray[n] = (i_temp + 1);
			}
		}

		if (i_bamRep > 500) { // THIS CALCULATES INFLECTION POINT OF REPLICATES SEPERATELY THEN AVERAGES.
			//drcNormScriptGeneratorString(i_loopCounter, i_bamFiles, i_truncDiff, turnoverTimes, normReadsVector, ($objName).b_realPeak, i_minIndex, i_bamRep, i_peakNumber);					
		} else {											
			drcVersitile(i_drcPar, i_maxThreads, i_loopCounter, i_bamFiles, turnoverTimes, ($objName).turnoverArr, ($objName).b_typeArr, i_minIndex, i_maxIndex, d_satthresh, b_forwardArr, b_reverseArr, b_peak1Arr, b_peak2Arr, b_dip1Arr, b_dip2Arr, b_model);
		}												

		// Progress
		std::cout << "[";
		int i_pos = i_barWidth * f_progress;
		for (int i = 0; i < i_barWidth; i++) {
			if (i <= i_pos) std::cout << "=";
			else std::cout << " ";
		}
		if (i == i_peakNumber) { // 100% complete
			std::cout << "] " << int(100) << " %\r";
			std::cout.flush();
		} else {
			std::cout << "] " << int((f_progress+0.01) * 100.0) << " %\r";
			std::cout.flush();
		}

		f_progress = d_progiterator/i_peakNumber;
		d_progiterator++;
		i_loopCounter++;
	} // end of object loop
	std::cout << "\nAnalysing results." << std::endl;
	
	inflectionWriterVersitile(i_maxThreads, i_peakNumber, s_drcFI_vec, s_drcFH_vec, s_drcFU_vec, s_drcFL_vec, s_drcFE_vec, s_drcRI_vec, s_drcRH_vec, s_drcRU_vec, s_drcRL_vec, s_drcRE_vec, s_type_vec);


	if (b_model == false) { // fix the s_type_vec values from undefined to rises/fall
		for (int i = 0; i < i_peakNumber; i++) {
			if (s_drcFH_vec[i] == "-") { s_type_vec[i] = "fall"; }
			else { s_type_vec[i] = "rise"; }
		}
	}

	/** 
	THERE IS A PROBLEM WITH THE MODELLING ALGORITHM. SOME LOCI ARE MODELLED AS A CERTAIN TYPE BUT FAIL TO BEHAVE THAT WAY.
	FOR EXAMPLE, A LOCI CAN BE MODELLED AS A RISE BUT HAVE A POSITIVE HILL'S COEFFICIENT. OR VALLEYS/HILLS COULD HAVE 
	2 HILLS COEFFICIENT OF THE SAME SIGN. EASY SOLUTION WOULD BE TO SIMPLY CHECK DATA ARRAY AFTER MODELLING IS COMPLETE AND 
	WRITE OVER THESE LOCI WITH ELIMINATED
	**/
	//if (b_model) {
		// ELIMINATE BAD LOCI
		for (int i = 0; i < i_peakNumber; i++) {
			if ( (s_drcFH_vec[i] == "NaN") || (s_drcRH_vec[i] == "NaN") ) {
				s_type_vec[i] = "eliminated";
				s_drcFI_vec[i] = "NaN"; 
				s_drcFH_vec[i] = "NaN";
				s_drcFU_vec[i] = "NaN"; 
				s_drcFL_vec[i] = "NaN";
				s_drcFE_vec[i] = "NaN";		
				s_drcRI_vec[i] = "NaN"; 
				s_drcRH_vec[i] = "NaN"; 
				s_drcRU_vec[i] = "NaN"; 
				s_drcRL_vec[i] = "NaN"; 
				s_drcRE_vec[i] = "NaN"; 
			} else if ( ((s_type_vec[i] == "rise") && (stod(s_drcFH_vec[i]) > 0)) || // bad rise
			     ((s_type_vec[i] == "fall") && (stod(s_drcRH_vec[i]) < 0)) || // bad fall
			     ((s_type_vec[i] == "hill") && ((stod(s_drcRH_vec[i]) < 0) || (stod(s_drcFH_vec[i]) > 0))) ||     // bad hill
			     ((s_type_vec[i] == "valley") && ((stod(s_drcRH_vec[i]) < 0) || (stod(s_drcFH_vec[i]) > 0))) ) {  // bad valley
					s_type_vec[i] = "eliminated";
					s_drcFI_vec[i] = "NaN"; 
					s_drcFH_vec[i] = "NaN";
					s_drcFU_vec[i] = "NaN"; 
					s_drcFL_vec[i] = "NaN";	
					s_drcFE_vec[i] = "NaN";		
					s_drcRI_vec[i] = "NaN"; 
					s_drcRH_vec[i] = "NaN"; 
					s_drcRU_vec[i] = "NaN"; 
					s_drcRL_vec[i] = "NaN"; 
					s_drcRE_vec[i] = "NaN"; 
			} else {}
		}
	//}

	// USE F VALUE TO ADJUST INFLECTION POINT 
	for (int i = 0; i < i_peakNumber; i++) {
		if (b_model) {
			if (s_type_vec[i] == "rise") { 
				double d = stod(s_drcFU_vec[i]);
				double c;
				if (s_drcFL_vec[i] == "forced to 0") { c = 0; } else { c = stod(s_drcFL_vec[i]); }
				double b = stod(s_drcFH_vec[i]);
				double e = stod(s_drcFI_vec[i]);
				double f = stod(s_drcFE_vec[i]);
				double y = (d-c)/2;
				double d_adjustedIP = ((log(pow(((d-c)/(y-c)),(1/f))-1))/b)+e;
				if ( isnan(d_adjustedIP) || isinf(d_adjustedIP) ) { }
				else { s_drcFI_vec[i] = std::to_string(d_adjustedIP); }
			} else if (s_type_vec[i] == "fall") { 
				double d = stod(s_drcRU_vec[i]);
				double c;
				if (s_drcRL_vec[i] == "forced to 0") { c = 0; } else { c = stod(s_drcRL_vec[i]); }
				double b = stod(s_drcRH_vec[i]);
				double e = stod(s_drcRI_vec[i]);
				double f = stod(s_drcRE_vec[i]);
				double y = (d-c)/2;
				double d_adjustedIP = ((log(pow(((d-c)/(y-c)),(1/f))-1))/b)+e;
				if ( isnan(d_adjustedIP) || isinf(d_adjustedIP) ) { }
				else { s_drcRI_vec[i] = std::to_string(d_adjustedIP); }
			} else if (s_type_vec[i] == "hill") { 
				double d = stod(s_drcFU_vec[i]);
				double c;
				if (s_drcFL_vec[i] == "forced to 0") { c = 0; } else { c = stod(s_drcFL_vec[i]); }
				double b = stod(s_drcFH_vec[i]);
				double e = stod(s_drcFI_vec[i]);
				double f = stod(s_drcFE_vec[i]);
				double y = (d-c)/2;
				double d_adjustedIP = ((log(pow(((d-c)/(y-c)),(1/f))-1))/b)+e;
				if ( isnan(d_adjustedIP) || isinf(d_adjustedIP) ) { }
				else { s_drcFI_vec[i] = std::to_string(d_adjustedIP); }

				d = stod(s_drcRU_vec[i]);
				if (s_drcRL_vec[i] == "forced to 0") { c = 0; } else { c = stod(s_drcRL_vec[i]); }
				b = stod(s_drcRH_vec[i]);
				e = stod(s_drcRI_vec[i]);
				f = stod(s_drcRE_vec[i]);
				y = (d-c)/2;
				d_adjustedIP = ((log(pow(((d-c)/(y-c)),(1/f))-1))/b)+e;
				if ( isnan(d_adjustedIP) || isinf(d_adjustedIP) ) { }
				else { s_drcRI_vec[i] = std::to_string(d_adjustedIP); }

			} else if (s_type_vec[i] == "valley") { 
				double d = stod(s_drcFU_vec[i]);
				double c;
				if (s_drcFL_vec[i] == "forced to 0") { c = 0; } else { c = stod(s_drcFL_vec[i]); }
				double b = stod(s_drcFH_vec[i]);
				double e = stod(s_drcFI_vec[i]);
				double f = stod(s_drcFE_vec[i]);
				double y = (d-c)/2;
				double d_adjustedIP = ((log(pow(((d-c)/(y-c)),(1/f))-1))/b)+e;
				if ( isnan(d_adjustedIP) || isinf(d_adjustedIP) ) { }
				else { s_drcFI_vec[i] = std::to_string(d_adjustedIP); }

				d = stod(s_drcRU_vec[i]);
				if (s_drcRL_vec[i] == "forced to 0") { c = 0; } else { c = stod(s_drcRL_vec[i]); }
				b = stod(s_drcRH_vec[i]);
				e = stod(s_drcRI_vec[i]);
				f = stod(s_drcRE_vec[i]);
				y = (d-c)/2;
				d_adjustedIP = ((log(pow(((d-c)/(y-c)),(1/f))-1))/b)+e;
				if ( isnan(d_adjustedIP) || isinf(d_adjustedIP) ) { }
				else { s_drcRI_vec[i] = std::to_string(d_adjustedIP); }

			} else if ( (s_type_vec[i] == "undefined") && (s_drcRH_vec[i] == "-") ) { // undefined rise
				double d = stod(s_drcFU_vec[i]);
				double c;
				if (s_drcFL_vec[i] == "forced to 0") { c = 0; } else { c = stod(s_drcFL_vec[i]); }
				double b = stod(s_drcFH_vec[i]);
				double e = stod(s_drcFI_vec[i]);
				double f = stod(s_drcFE_vec[i]);
				double y = (d-c)/2;
				double d_adjustedIP = ((log(pow(((d-c)/(y-c)),(1/f))-1))/b)+e;
				if ( isnan(d_adjustedIP) || isinf(d_adjustedIP) ) { }
				else { s_drcFI_vec[i] = std::to_string(d_adjustedIP); }
			} else if ( (s_type_vec[i] == "undefined") && (s_drcFH_vec[i] == "-") ) { // undefined fall
				double d = stod(s_drcRU_vec[i]);
				double c;
				if (s_drcRL_vec[i] == "forced to 0") { c = 0; } else { c = stod(s_drcRL_vec[i]); }
				double b = stod(s_drcRH_vec[i]);
				double e = stod(s_drcRI_vec[i]);
				double f = stod(s_drcRE_vec[i]);
				double y = (d-c)/2;
				double d_adjustedIP = ((log(pow(((d-c)/(y-c)),(1/f))-1))/b)+e;
				if ( isnan(d_adjustedIP) || isinf(d_adjustedIP) ) { }
				else { s_drcRI_vec[i] = std::to_string(d_adjustedIP); }
			} else {}
		} else { // nomodel
			if (s_drcFH_vec[i] == "-") { // reverse peak (fall)
				double d = stod(s_drcRU_vec[i]);
				double c;
				if (s_drcRL_vec[i] == "forced to 0") { c = 0; } else { c = stod(s_drcRL_vec[i]); }
				double b = stod(s_drcRH_vec[i]);
				double e = stod(s_drcRI_vec[i]);
				double f = stod(s_drcRE_vec[i]);
				double y = (d-c)/2;
				double d_adjustedIP = ((log(pow(((d-c)/(y-c)),(1/f))-1))/b)+e;
				if ( isnan(d_adjustedIP) || isinf(d_adjustedIP) ) { }
				else { s_drcRI_vec[i] = std::to_string(d_adjustedIP); }
			} else { // forward peak (rise)
				double d = stod(s_drcFU_vec[i]);
				double c;
				if (s_drcFL_vec[i] == "forced to 0") { c = 0; } else { c = stod(s_drcFL_vec[i]); }
				double b = stod(s_drcFH_vec[i]);
				double e = stod(s_drcFI_vec[i]);
				double f = stod(s_drcFE_vec[i]);
				double y = (d-c)/2;
				double d_adjustedIP = ((log(pow(((d-c)/(y-c)),(1/f))-1))/b)+e;
				if ( isnan(d_adjustedIP) || isinf(d_adjustedIP) ) { }
				else { s_drcFI_vec[i] = std::to_string(d_adjustedIP); }
			}
		}
	}


	// TRANSFER DATA TO dataArray
	for (int i = 0; i < i_peakNumber; i++) {
		dataArray[i+1][2*i_bamFiles + 3]  = s_type_vec[i]; 
		dataArray[i+1][2*i_bamFiles + 4]  = s_drcFI_vec[i];  // F as in forward
		dataArray[i+1][2*i_bamFiles + 5]  = s_drcFH_vec[i];
		dataArray[i+1][2*i_bamFiles + 6]  = s_drcFU_vec[i]; 
		dataArray[i+1][2*i_bamFiles + 7]  = s_drcFL_vec[i];	
		dataArray[i+1][2*i_bamFiles + 8]  = s_drcRI_vec[i];  // R as in reverse
		dataArray[i+1][2*i_bamFiles + 9]  = s_drcRH_vec[i]; 
		dataArray[i+1][2*i_bamFiles + 10] = s_drcRU_vec[i]; 
		dataArray[i+1][2*i_bamFiles + 11] = s_drcRL_vec[i]; 
	}

	// GRAPHING COMPONENT
	std::cout << "Writing R script." << std::endl;

	int i_undef      = 0;
	int i_elim       = 0;
	int i_valley     = 0;
	int i_hill       = 0;
	int i_fall       = 0;
	int i_rise       = 0;
	int i_undefRise  = 0;
	int i_undefFall  = 0;
	int i_increase   = 0;
	int i_decrease   = 0;

	for (int i = 0; i < i_peakNumber; i++) {
		if (dataArray[i+1][2*i_bamFiles + 3]   == "eliminated") { i_elim++;   }
		if (dataArray[i+1][2*i_bamFiles + 3]   == "undefined")  { i_undef++;  }
		if (dataArray[i+1][2*i_bamFiles + 3]   == "valley")     { i_valley++; }
		if (dataArray[i+1][2*i_bamFiles + 3]   == "hill")       { i_hill++;   }
		if (dataArray[i+1][2*i_bamFiles + 3]   == "fall")       { i_fall++;   }
		if (dataArray[i+1][2*i_bamFiles + 3]   == "rise")       { i_rise++;   }
		if ( (dataArray[i+1][2*i_bamFiles + 3] == "undefined") && (dataArray[i+1][2*i_bamFiles + 5]) != "-" )  { i_undefRise++;  }
		if ( (dataArray[i+1][2*i_bamFiles + 3] == "undefined") && (dataArray[i+1][2*i_bamFiles + 5]) == "-" )  { i_undefFall++;  }
	}
	i_increase = i_undefRise + i_rise + i_valley + i_hill;
	i_decrease = i_undefFall + i_fall + i_valley + i_hill;

/**
	XXX ============================= XXX
	XXX FIX GRAPHING COMPONENTS BELOW XXX =========================================================== XXX
	XXX ============================= XXX
	XXX TODO:

	FIX GRAPHS:	
		1. Giemsa staining and scales on ideograms
	CORE PROGRAM CHANGES:
		1. change last part of stdev.cpp to import data from file 
		2. make a function that checks if all dependencies are available 
		(expand pkgCheck.cpp and create somthing before to check if R, samtools and bedtools is available - make a external software tester?)

		find centromeres for sacCer3 ideogram 
		fix hg19 ideogram - probably hg38 also... (double check centromere positions)

	XXX ============================================================================================== XXX
**/

	// pie plot XXX FOR nomodel, FIX SO THAT ONLY FORWARD AND REVERSE ARE DISPLAYED XXX
	rPie(s_rPltsName, i_peakNumber, i_undef, i_elim, i_valley, i_hill, i_fall, i_rise, b_model);

	// min/max bar chart
	rBar(s_rPltsName, i_bamFiles, turnoverTimes, maxIndexArray, minIndexArray, i_peakNumber);

	// normalized heamap of all peaks
	rHeatmap(s_rPltsName, i_peakNumber, i_bamFiles, dataArray, turnoverTimes);

	// average profiles
	rProfiles(s_rPltsName, i_peakNumber, i_bamFiles, dataArray, b_model, d_satthresh, turnoverTimes, i_valley, i_hill, i_fall, i_rise, i_undefRise, i_undefFall);

	// hill / valley quality charts
	if (b_model) {
		// scatterplot boxplot of log2(rise/fall) of upper aymptotes in hills and valleys.
		// boxplot / bar chart of hills coefficient for rises, hill rises, valley rises, and undef rises.
		rUpperA(s_rPltsName, i_peakNumber, i_bamFiles, dataArray, i_valley, i_hill);
	}

	// boxplot of gene features and ideogram heat map 		
	if ( b_genomeFlagCall == true ) {
		if (i_increase > 0) {
			rBoxplot(turnoverTimes, s_rPltsName,i_peakNumber,i_genomeFileCount,i_bamFiles,dataArray,std::string("forward"));
		} 
		if (i_decrease > 0) {
			rBoxplot(turnoverTimes, s_rPltsName,i_peakNumber,i_genomeFileCount,i_bamFiles,dataArray,std::string("reverse"));
		} 
		if ( s_genome == "hg19") {
			if (i_increase > 0) {
				rIdeohg19(s_rPltsName,i_bamFiles,i_peakNumber,dataArray, std::string("forward"));
			} 
			if (i_decrease > 0) {
				rIdeohg19(s_rPltsName,i_bamFiles,i_peakNumber,dataArray, std::string("reverse")); 
			}
		}		
		if ( s_genome == "hg38") {
			if (i_increase > 0) {
				rIdeohg38(s_rPltsName,i_bamFiles,i_peakNumber,dataArray, std::string("forward"));
			} 
			if (i_decrease > 0) {
				rIdeohg38(s_rPltsName,i_bamFiles,i_peakNumber,dataArray, std::string("reverse")); 
			} 
		}		
		if ( s_genome == "mm9") {
			if (i_increase > 0) {
				rIdeomm9(s_rPltsName,i_bamFiles,i_peakNumber,dataArray, std::string("forward"));
			} 
			if (i_decrease > 0) {
				rIdeomm9(s_rPltsName,i_bamFiles,i_peakNumber,dataArray, std::string("reverse")); 
			}
		}
		if ( s_genome == "mm10") {
			if (i_increase > 0) {
				rIdeomm10(s_rPltsName,i_bamFiles,i_peakNumber,dataArray, std::string("forward"));
			} 
			if (i_decrease > 0) {
				rIdeomm10(s_rPltsName,i_bamFiles,i_peakNumber,dataArray, std::string("reverse")); 
			} 
		}
		if ( s_genome == "dm6") {
			if (i_increase > 0) {
				rIdeodm6(s_rPltsName,i_bamFiles,i_peakNumber,dataArray, std::string("forward"));
			} 
			if (i_decrease > 0) {
				rIdeodm6(s_rPltsName,i_bamFiles,i_peakNumber,dataArray, std::string("reverse")); 
			} 
		}
		if ( s_genome == "dm3") {
			if (i_increase > 0) {
				rIdeodm3(s_rPltsName,i_bamFiles,i_peakNumber,dataArray, std::string("forward"));
			} 
			if (i_decrease > 0) {
				rIdeodm3(s_rPltsName,i_bamFiles,i_peakNumber,dataArray, std::string("reverse")); 
			} 
		}	
		if ( s_genome == "ce11") {
			if (i_increase > 0) {
				rIdeoce11(s_rPltsName,i_bamFiles,i_peakNumber,dataArray, std::string("forward"));
			} 
			if (i_decrease > 0) {
				rIdeoce11(s_rPltsName,i_bamFiles,i_peakNumber,dataArray, std::string("reverse")); 
			}
		}
		if ( s_genome == "sacCer3") {
			if (i_increase > 0) {
				rIdeosacCer3(s_rPltsName,i_bamFiles,i_peakNumber,dataArray, std::string("forward"));
			} 
			if (i_decrease > 0) {
				rIdeosacCer3(s_rPltsName,i_bamFiles,i_peakNumber,dataArray, std::string("reverse"));
			} 
		}
	}

	// stdev graphs
	if (i_bamRep > 1) 
		rStdev(s_rPltsName, i_bamFiles, i_peakNumber, dataArray, turnoverTimes, i_undefRise, i_undefFall, i_valley, i_hill, i_fall, i_rise); 

	// Scatter bars
	if (i_increase > 0) 
		rScat(s_rPltsName, i_bamFiles, turnoverTimes, i_peakNumber, dataArray, i_undefRise, i_undefFall, i_valley, i_hill, i_fall, i_rise, b_model, std::string("forward"));
	if (i_decrease > 0) 
		rScat(s_rPltsName, i_bamFiles, turnoverTimes, i_peakNumber, dataArray, i_undefRise, i_undefFall, i_valley, i_hill, i_fall, i_rise, b_model, std::string("reverse"));


	// 3D GRAPHS
	if (b_genelistFlagCall) {
		// R script file
		std::string s_rGenePltsName = s_name + ".tdca3Dgenes.R";
		std::vector<double> shrunkGeneVector;

		// 3d scatter plots
		if (i_bamRep > 1 && i_bamFiles < 6) {
			r3Dadjusted(s_genelistFile, s_rGenePltsName, i_bamFiles, turnoverTimes, i_validgenes, normGeneAveVector, validGeneVector);
		} else if (i_bamRep == 1 && i_bamFiles < 6) {
			r3Dadjusted(s_genelistFile, s_rGenePltsName, i_bamFiles, turnoverTimes, i_validgenes, normGeneVector, validGeneVector);
		} else if (i_bamRep > 1 && i_bamFiles > 5) {
			geneVectorShrinkage(i_validgenes, i_bamFiles, normGeneAveVector, shrunkGeneVector);
			r3Dadjusted(s_genelistFile, s_rGenePltsName, i_bamFiles, turnoverTimes, i_validgenes, shrunkGeneVector, validGeneVector);
		} else if (i_bamRep == 1 && i_bamFiles > 5) {
			geneVectorShrinkage(i_validgenes, i_bamFiles, normGeneVector, shrunkGeneVector);
			r3Dadjusted(s_genelistFile, s_rGenePltsName, i_bamFiles, turnoverTimes, i_validgenes, shrunkGeneVector, validGeneVector);
		} else {
			std::cout << "Error with 3d gene graphing." << std::endl;	
		}

		// run r script
		std::string s_rGeneScript = "Rscript " + s_rGenePltsName; // temp file
		try {
			std::system(s_rGeneScript.c_str());
		}
		catch(...){
			std::cerr << "Cannot access R" << std::endl;
			std::exit(0);
		}
		// change name of Rplot.pdf
		std::string s_changeName = "mv Rplots.pdf " + s_name + ".tdca3Dgenes.pdf";
		try {
			std::system(s_changeName.c_str());
		}
		catch(...){
			std::cerr << "Cannot change name of Rplot" << std::endl;
			std::exit(0);
		}
	}
	rEnd(s_rPltsName, s_rPltsPDF, b_genomeFlagCall, i_bamRep, b_model, i_increase, i_decrease); 


	std::cout << "Printing data to " << s_name << ".tdca.txt" << std::endl;
	// rows = i_peakNumber + 1; columns	= 2*i_bamFiles + i_genomeFileCount + 11
	std::ofstream outputFile;
	outputFile.open(s_name + ".tdca.txt"); 

	for (int i = 0; i < i_peakNumber; i++) {
		for (int j = 0; j < (i_dataArrCols-1); j++) 
			outputFile << dataArray[i][j] << '\t';

		outputFile << dataArray[i][(i_dataArrCols-1)] << '\n';		
	}

	for (int i = 0; i < (i_dataArrCols-1); i++) 
		outputFile << dataArray[i_peakNumber][i] << '\t';
	outputFile << dataArray[i_peakNumber][(i_dataArrCols-1)];


	validGeneVector.clear();
	normGeneAveVector.clear();
	normGeneVector.clear();
	dataArray.clear();
	outputFile.close();

} //main
