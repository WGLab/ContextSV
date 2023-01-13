//
// Created by jperdomo on 1/8/2023.
//

#include "cli.h"
#include "bam_reader.h"

#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>
#include <stdexcept>

cli::cli() = default;

// Check if the filepath exists
bool cli::fileExists(const std::string &name) {
	if (FILE *file = fopen(name.c_str(), "r")) {
		fclose(file);
		return true;
	} else {
		return false;
	}
}

// Get the command argument input
std::string cli::getCmdOption(char** begin, char** end, const std::string& option)
{
	std::string result("");
	char ** itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
	{
		result = *itr;
		//return *itr;
	}
	return result;
}

// Check if the input argument is provided
bool cli::cmdOptionExists(char** begin, char** end, const std::string& option)
{
	return std::find(begin, end, option) != end;
}

std::string cli::getInputFilepath() {
    return input_filepath;
}

// Parse input arguments
int cli::parse(int argc, char **argv) {
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
		// Get the input BAM file
		std::string filename = getCmdOption(argv, argv + argc, "--bam");

		if (fileExists(filename)) {
			this->input_filepath = filename;
			std::cout << "Input BAM = " << this->input_filepath << std::endl;
			
			// Set the success code
			exit_code = 0;
		} else {
			std::string err_str = "Input BAM does not exist: " + filename;
			throw std::invalid_argument(err_str);
		}
	}

	return exit_code;
}

// Run the CLI with input arguments
void cli::run()
{
	// Read the file
	std::string filepath = getInputFilepath();
	samFile * in_bam;
	bam_hdr_t * hdr;
	in_bam = sam_open(filepath.c_str(), "rb");
	hdr = sam_hdr_read(in_bam);
	std::cout << "SAMTools succeeded" << std::endl;
}

void cli::printHelpText() {
    int width = 20;
    std::ios_base::fmtflags flags = std::cout.flags();

    std::cout << std::left
              << std::setw(width) << "\nPositional arguments:\n"
              << std::setw(width) << "-i, --input" << std::setw(20) << "BAM file input"
              << std::endl
              << std::setw(width) << "-o, --output" << std::setw(20) << "Output file directory"
              << std::endl
              << std::setw(width) << "\nOptional arguments:\n"
              << std::setw(width) << "-h, --help" << std::setw(20) << "Show this help message and exit\n"
              << std::endl;

    std::cout.flags(flags);
}
