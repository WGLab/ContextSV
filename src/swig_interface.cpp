//
// Created by jperdomo on 1/8/2023.
//

#include "swig_interface.h"
#include "contextsv.h"

/// @cond
#include <iostream>
/// @endcond


// Run the CLI with the given parameters
int run(std::string bam_fp, std::string ref_fp, std::string snps_fp, std::string outdir, std::string region, std::string chr_cov, std::string pfb_fp)
{
	// Create the input_data parameters
	InputData input_data;
	input_data.setBAMFilepath(bam_fp);
	input_data.setRefGenome(ref_fp);
	input_data.setSNPFilepath(snps_fp);
	input_data.setOutputDir(outdir);
	input_data.setRegion(region);
	input_data.setChrCov(chr_cov);
	input_data.setPFBFilepath(pfb_fp);

	// Run ContextSV
	ContextSV caller_obj(input_data);
	try
	{	
		caller_obj.run();
	}

	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
		return -1;
	}

	return 0;
}
