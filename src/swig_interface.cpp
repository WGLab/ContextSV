//
// Created by jperdomo on 1/8/2023.
//

#include "swig_interface.h"
#include "contextsv.h"

#include <iostream>


// Run the CLI with the given parameters
int run(std::string bam, std::string snps, std::string outdir, std::string region)
{
	// Create the common parameters
	Common common;
	common.setBAMFilepath(bam);
	common.setSNPFilepath(snps);
	common.setOutputDir(outdir);
	common.setRegion(region);

	// Run ContextSV
	ContextSV caller_obj(common);
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
