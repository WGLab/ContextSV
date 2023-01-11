//
// Created by jperdomo on 1/8/2023.
//

#include "cli.h"
#include "bam_reader.h"

#include <htslib/sam.h>
#include <algorithm>
#include <string>
#include <iostream>
#include <iomanip>

cli::cli() = default;

// Check if a file exists
bool cli::fileExists(const std::string &name) {
	if (FILE *file = fopen(name.c_str(), "r")) {
		fclose(file);
		return true;
	} else {
		return false;
	}
}

// Get the argument option
char* cli::getCmdOption(char** begin, char** end, const std::string& option)
{
	char ** itr = std::find(begin, end, option);
	if (itr != end && ++itr != end)
	{
		return *itr;
	}
	return nullptr;
}

// Check if the option was provided
bool cli::cmdOptionExists(char** begin, char** end, const std::string& option)
{
	return std::find(begin, end, option) != end;
}

// Get input filepath
std::string cli::getInputFilepath() {
    return input_filepath;
}

// Parse arguments
void cli::parse(int argc, char **argv) {

    // Print help text
	if (argc == 1 || cmdOptionExists(argv, argv+argc, "-h") || cmdOptionExists(argv, argv+argc, "--help"))
	{
		// No arguments specified
		printHelpText();
	} else
	{
		// Get the input BAM file
		std::string filename = getCmdOption(argv, argv + argc, "--bam");
		if (fileExists(filename)) {
			this->input_filepath = filename;
			std::cout << "Input BAM = " << this->input_filepath << std::endl;

			// Read the file
			samFile * in_bam;
			bam_hdr_t * hdr;
			in_bam = sam_open(filename.c_str(), "rb");
			hdr = sam_hdr_read(in_bam);
			std::cout << "SAMTools succeeded" << std::endl;
		} else {
			std::cerr << "Input BAM does not exist: " << filename << std::endl;
		}
	}
}

void cli::read()
{
}

// Print help text
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
