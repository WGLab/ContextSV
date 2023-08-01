//
// Created by jperdomo on 1/8/2023.
//

#include "cli.h"
#include "integrative_caller.h"

#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <htslib/sam.h>


// Run the CLI with the given parameters
int run(std::string bam, std::string snps, std::string outdir, std::string region)
{
	// Create the common parameters
	Common common;
	common.setBAMFilepath(bam);
	common.setSNPFilepath(snps);
	common.setOutputDir(outdir);
	common.setRegion(region);

	// Run the integrative caller
	IntegrativeCaller caller_obj(common);
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
