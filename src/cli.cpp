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


// Check if the filepath exists
bool CLI::fileExists(const std::string &name) {
	if (FILE *file = fopen(name.c_str(), "r")) {
		fclose(file);
		return true;
	} else {
		return false;
	}
}

// Get the command argument input
std::string CLI::getCmdOption(char** begin, char** end, const std::string& option)
{
	std::string result("");
	char ** itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
	{
		result = *itr;
	}
	return result;
}

// Check if the input argument is provided
bool CLI::cmdOptionExists(char** begin, char** end, const std::string& option)
{
	return std::find(begin, end, option) != end;
}

// Parse input arguments
int CLI::parse(int argc, char **argv) {
	// Initialize the common parameters
	this->common = Common();

	// Exit code 0 indicates parameters were successfully set
	int exit_code = 0;

    // Print help text
	if (argc == 1 || cmdOptionExists(argv, argv+argc, "-h") || cmdOptionExists(argv, argv+argc, "--help"))
	{
		// No arguments specified
		printHelpText();
	}
	
	else
	{
		// Get the output directory
		std::string output_dir = getCmdOption(argv, argv + argc, "-o");
		if (output_dir.empty()) {
			output_dir = getCmdOption(argv, argv + argc, "--out");
		}
		if (!output_dir.empty()) {
			this->common.set_output_dir(output_dir);
			std::cout << "Output directory = " << output_dir << std::endl;

		} else {
			std::string err_str = "Output directory not specified.";
			throw std::invalid_argument(err_str);
		}

		// Get the reference genome file
		std::string ref_filename = getCmdOption(argv, argv + argc, "--ref");
		if (fileExists(ref_filename)) {
			common.set_ref_filepath(ref_filename);
			std::cout << "Reference genome file = " << ref_filename << std::endl;

		} else {
			std::string err_str = "File " + ref_filename + " does not exist.";
			throw std::invalid_argument(err_str);
		}


		// Get the bam file
		std::string bam_filename = getCmdOption(argv, argv + argc, "--bam");
		if (fileExists(bam_filename)) {
			common.set_bam_filepath(bam_filename);
			std::cout << "Alignment file = " << bam_filename << std::endl;
			
		} else {
			std::string err_str = "BAM file does not exist: " + bam_filename;
			throw std::invalid_argument(err_str);
		}

		// Get the region to analyze
		std::string region = getCmdOption(argv, argv + argc, "--region");
		if (!region.empty()) {
			//this->region = region;
			common.set_region(region);
			std::cout << "Region = " << region << std::endl;
		}

		// Get the window size
		std::string window_size = getCmdOption(argv, argv + argc, "--window-size");
		if (!window_size.empty()) {
			common.set_window_size(std::stoi(window_size));
			std::cout << "Window size = " << window_size << std::endl;
		}

		// Get the SNP VCF file
		std::string snp_vcf_filename = getCmdOption(argv, argv + argc, "--snp-vcf");
		if (fileExists(snp_vcf_filename)) {
			common.set_snp_vcf_filepath(snp_vcf_filename);
			std::cout << "SNP VCF file = " << snp_vcf_filename << std::endl;

		} else if (!snp_vcf_filename.empty()) {
			std::string err_str = "SNP VCF file does not exist: " + snp_vcf_filename;
			throw std::invalid_argument(err_str);
		}
	}

	return exit_code;
}

// Run the CLI with input arguments
int CLI::run()
{
	// Run the integrative caller
	IntegrativeCaller caller_obj(this->common);
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

void CLI::printHelpText() {
    int width = 20;
    std::ios_base::fmtflags flags = std::cout.flags();

    std::cout << std::left
              << std::setw(width) << "\nPositional arguments:\n"
			  << std::setw(width) << "--ref" << std::setw(20) << "reference genome fasta file"
              << std::setw(width) << "--bam" << std::setw(20) << "alignment file in BAM format"
              << std::endl
              << std::setw(width) << "-o, --output" << std::setw(20) << "output file directory"
              << std::endl
			  << std::setw(width) << "--region" << std::setw(20) << "region to analyze"
              << std::setw(width) << "\nOptional arguments:\n"
			  << std::setw(width) << "--window-size" << std::setw(20) << "window size (default = 10000)"
			  << std::endl
              << std::setw(width) << "-h, --help" << std::setw(20) << "Show this help message and exit\n"
              << std::endl;

    std::cout.flags(flags);
}
