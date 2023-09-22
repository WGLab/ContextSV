#include "fasta_query.h"

#include <string.h>
#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <sstream>


int FASTAQuery::setFilepath(std::string fasta_filepath)
{
    if (fasta_filepath == "")
    {
        std::cout << "No FASTA filepath provided" << std::endl;
        return 1;
    }

    this->fasta_filepath = fasta_filepath;

    std::cout << "Reading FASTA file " << fasta_filepath << std::endl;

    // Parse the FASTA file
    std::ifstream fasta_file(fasta_filepath);
    if (!fasta_file.is_open())
    {
        std::cout << "Could not open FASTA file " << fasta_filepath << std::endl;
        exit(1);
    }

    // Get the sequence
    std::map<std::string, std::string> chr_to_seq;
    std::string current_chr = "";
    std::string sequence = "";
    std::string line_str = "";

    while (std::getline(fasta_file, line_str))
    {
        // Check if the line is a header
        if (line_str[0] == '>')
        {
            // Header line, indicating a new chromosome
            // Store the previous chromosome and sequence
            if (current_chr != "")
            {
                chr_to_seq[current_chr] = sequence;
                std::cout << "Read chromosome " << current_chr << " with length " << sequence.length() << std::endl;
                sequence = "";  // Reset the sequence
            }

            // Get the new chromosome
            current_chr = line_str.substr(1);

            // Remove the description
            size_t space_pos = current_chr.find(" ");
            if (space_pos != std::string::npos)
            {
                current_chr.erase(space_pos);
            }

            // Check if the chromosome is already in the map
            if (chr_to_seq.find(current_chr) != chr_to_seq.end())
            {
                std::cout << "Duplicate chromosome " << current_chr << std::endl;
                exit(1);
            }
        } else {
            // Sequence line
            sequence += line_str;
        }
    }

    // Add the last chromosome at the end of the file
    if (current_chr != "")
    {
        chr_to_seq[current_chr] = sequence;
    }

    // Close the file
    std::cout << "Closing FASTA file..." << std::endl;
    fasta_file.close();

    // Set the map
    this->chr_to_seq = chr_to_seq;

    std::cout << "Done." << std::endl;

    return 0;
}

std::string FASTAQuery::getFilepath()
{
    return this->fasta_filepath;
}

// Function to get the reference sequence at a given position range
std::string FASTAQuery::query(std::string chr, int pos_start, int pos_end)
{
    std::cout << "Querying " << chr << ":" << pos_start << "-" << pos_end << std::endl;

    // Check if a FASTA file has been set
    if (this->fasta_filepath == "")
    {
        std::cout << "No FASTA file set" << std::endl;
        return "";
    }

    // Check if the chromosome is in the map
    if (this->chr_to_seq.find(chr) == this->chr_to_seq.end())
    {
        std::cout << "Chromosome " << chr << " not found in FASTA file" << std::endl;
        return "";
    }
    
    //std::cout << "Querying " << chr << ":" << pos_start << "-" << pos_end << std::endl;

    // Get the sequence
    std::string sequence = this->chr_to_seq[chr];

    //std::cout << "Sequence length: " << sequence.length() << std::endl;

    // Get the substring
    std::string subsequence = sequence.substr(pos_start, pos_end - pos_start);

    //std::cout << "Subsequence length: " << subsequence.length() << std::endl;

    // If the subsequence is empty, return VCF missing data character
    if (subsequence == "")
    {
        return ".";
    }

    return subsequence;
}
