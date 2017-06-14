#include <iostream>    // input/output
#include <string>
#include <vector>

void normGeneVectorAverage(int i_bamFiles, int i_bamRep, int i_validgenes, std::vector<double> &normGeneVector, std::vector<double> &normGeneAveVector, double d_bin)
{

	for (int i = 0; i < i_validgenes; i++) {     			 // for all genes
		for( int k = 0; k < i_bamFiles; k++) {	 		 // for each turnover time
			for (int j = 0; j < d_bin; j++) {      		 // for d_bin bins 
			double d_sumValue = 0;	
				for( int l = 0; l < i_bamRep; l++) {	 // for each bam replicate
					d_sumValue += normGeneVector[j+(d_bin*i_validgenes*l*i_bamFiles)+(d_bin*i*i_bamFiles)+(d_bin*k)];
				}
				double d_aveValue = d_sumValue/i_bamRep;
				normGeneAveVector.push_back(d_aveValue);
			}
		}
	}
}
