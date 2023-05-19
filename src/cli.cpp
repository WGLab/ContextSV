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

CLI::CLI() = default;

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

std::string CLI::get_bam_filepath() {
    return this->bam_filepath;
}

std::string CLI::get_ref_filepath()
{
    return this->ref_filepath;
}

// Parse input arguments
int CLI::parse(int argc, char **argv) {

	// Exit code 0 indicates parameters were successfully set
	int exit_code(1);

    // Print help text
	if (argc == 1 || cmdOptionExists(argv, argv+argc, "-h") || cmdOptionExists(argv, argv+argc, "--help"))
	{
		// No arguments specified
		printHelpText();
	}
	
	else
	{
		// Get the reference genome file
		std::string ref_filename = getCmdOption(argv, argv + argc, "--ref");
		if (fileExists(ref_filename)) {
			this->ref_filepath = ref_filename;
			std::cout << "Reference genome file = " << this->ref_filepath << std::endl;
		} else {
			std::string err_str = "File " + ref_filename + " does not exist.";
			throw std::invalid_argument(err_str);
		}


		// Get the bam file
		std::string bam_filename = getCmdOption(argv, argv + argc, "--bam");
		if (fileExists(bam_filename)) {
			this->bam_filepath = bam_filename;
			std::cout << "Alignment file = " << this->bam_filepath << std::endl;
			
			// Set the success code
			exit_code = 0;
		} else {
			std::string err_str = "File " + bam_filename + " does not exist.";
			throw std::invalid_argument(err_str);
		}
	}

	return exit_code;
}

// Run the CLI with input arguments
int CLI::run()
{
	// Get the input arguments
	std::string bam_filepath = this->get_bam_filepath();
	std::string ref_filepath = this->get_ref_filepath();
	IntegrativeCaller caller_obj;
	try
	{	
		caller_obj.set_bam_filepath(bam_filepath);
		caller_obj.set_ref_filepath(ref_filepath);
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
              << std::setw(width) << "\nOptional arguments:\n"
              << std::setw(width) << "-h, --help" << std::setw(20) << "Show this help message and exit\n"
              << std::endl;

    std::cout.flags(flags);
}
