//
// Created by jperdomo on 1/8/2023.
//

#include "swig_interface.h"
#include "contextsv.h"

#include <iostream>


// Run the CLI with the given parameters
int run(std::string bam_fp, std::string ref_fp, std::string snps_fp, std::string outdir, std::string region)
{
	// Create the input_data parameters
	InputData input_data;
	input_data.setBAMFilepath(bam_fp);
	input_data.setRefGenome(ref_fp);
	input_data.setSNPFilepath(snps_fp);
	input_data.setOutputDir(outdir);
	input_data.setRegion(region);

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
