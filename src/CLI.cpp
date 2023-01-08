//
// Created by jperdomo on 1/8/2023.
//

#include "CLI.h"

#include <string>
#include <cstring>
#include <iostream>
#include <iomanip>

CLI::CLI() = default;

std::string CLI::getInputFilepath() {
    return input_filepath;
}

void CLI::parse(int argc, char **argv) {
    // Print help text
    if (argc == 1) {
        // No arguments specified
        printHelpText();

    } else if (argc == 2) {
        // Option specified
        if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
            printHelpText();
        } else if (strcmp(argv[1], "-i") == 0 || strcmp(argv[1], "--input") == 0) {
            // Input filepath provided
            input_filepath = "test";
        } else {
            // Other arguments
        }
    }
}

void CLI::printHelpText() {
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

bool CLI::file_exists(const std::string &name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}
