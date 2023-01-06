//
// Created by perdomoj on 1/6/2023.
//

#include "ContextSV.h"

#include <iostream>
#include <iomanip>
#include <cstring>

// TODO: Move these to CLI file
void printHelpText() {
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

bool file_exists(const std::string& name) {
	if (FILE *file = fopen(name.c_str(), "r")) {
		fclose(file);
		return true;
	} else {
		return false;
	}
}

int main(int argc, char *argv[]){
    std::string bam_filepath;

    // Print help text
    if (argc == 1) {
        // No arguments specified
		printHelpText();

    } else if (argc == 2) {
		// Option specified
		if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
			printHelpText();
		} else {
			// Other arguments
		}
	}

    return 0;
}
