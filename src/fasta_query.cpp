#include "fasta_query.h"

/// @cond
#include <string.h>
#include <iostream>
#include <unordered_map>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
/// @endcond


#include "utils.h"


int ReferenceGenome::setFilepath(std::string fasta_filepath)
{
    if (fasta_filepath == "")
    {
        std::cout << "No FASTA filepath provided" << std::endl;
        return 1;
    }

    this->fasta_filepath = fasta_filepath;

    // Parse the FASTA file
    std::ifstream fasta_file(fasta_filepath);
    if (!fasta_file.is_open())
    {
        std::cout << "Could not open FASTA file " << fasta_filepath << std::endl;
        exit(1);
    }

    // Get the chromosomes and sequences
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
                this->chromosomes.push_back(current_chr);  // Add the chromosome to the list
                this->chr_to_seq[current_chr] = sequence;  // Add the sequence to the map
                this->chr_to_length[current_chr] = sequence.length();  // Add the sequence length to the map
                sequence = "";  // Reset the sequence
            }

            current_chr = line_str.substr(1);  // Remove the '>' character

            // Remove the description
            size_t space_pos = current_chr.find(" ");
            if (space_pos != std::string::npos)
            {
                current_chr.erase(space_pos);
            }
        } else {
            sequence += line_str;
        }
    }

    // Add the last chromosome at the end of the file
    if (current_chr != "")
    {
        this->chromosomes.push_back(current_chr);  // Add the chromosome to the list
        this->chr_to_seq[current_chr] = sequence;  // Add the sequence to the map
        this->chr_to_length[current_chr] = sequence.length();  // Add the sequence length to the map
    }

    fasta_file.close();
    std::sort(this->chromosomes.begin(), this->chromosomes.end());

    return 0;
}

std::string ReferenceGenome::getFilepath() const
{
    return this->fasta_filepath;
}

// Function to get the reference sequence at a given position range
std::string_view ReferenceGenome::query(const std::string& chr, uint32_t pos_start, uint32_t pos_end) const
{   
    // Convert positions from 1-indexed (reference) to 0-indexed (string indexing)
    pos_start--;
    pos_end--;

    // Ensure that the end position is not larger than the chromosome length
    const std::string& sequence = this->chr_to_seq.at(chr);
    if (pos_end >= sequence.length() || pos_start > pos_end)
    {
        return {};
    }

    return std::string_view(sequence).substr(pos_start, (pos_end - pos_start) + 1);
}

// Function to compare the reference sequence at a given position range
bool ReferenceGenome::compare(const std::string& chr, uint32_t pos_start, uint32_t pos_end, const std::string& compare_seq, float match_threshold) const
{    
    // Convert positions from 1-indexed (reference) to 0-indexed (string indexing)
    pos_start--;
    pos_end--;

    // Ensure that the end position is not larger than the chromosome length
    const std::string& sequence = this->chr_to_seq.at(chr);
    if (pos_end >= sequence.length() || pos_start >= pos_end)
    {
        return {};
    }

    std::string_view subseq = std::string_view(sequence).substr(pos_start, pos_end - pos_start + 1);
    if (subseq.length() != compare_seq.length())
    {
        printError("ERROR: Sequence lengths do not match for comparison");
        return false;
    }

    // Calculate the match rate
    size_t num_matches = 0;
    for (size_t i = 0; i < subseq.length(); i++)
    {
        if (subseq[i] == compare_seq[i])
        {
            num_matches++;
        }
    }
    float match_rate = (float)num_matches / (float)subseq.length();

    return match_rate >= match_threshold;
}

// Function to get the chromosome contig lengths in VCF header format
std::string ReferenceGenome::getContigHeader() const
{
    std::shared_lock<std::shared_mutex> lock(this->shared_mutex);
    std::string contig_header = "";

    // Sort the chromosomes
    std::vector<std::string> chromosomes;
    for (auto const& chr_seq : this->chr_to_seq)
    {
        chromosomes.push_back(chr_seq.first);
    }
    std::sort(chromosomes.begin(), chromosomes.end());
    for (auto const& chr : chromosomes)
    {
        // Add the contig header line
        contig_header += "##contig=<ID=" + chr + ",length=" + std::to_string(this->chr_to_seq.at(chr).length()) + ">\n";
    }

    // Remove the last newline character
    contig_header.pop_back();

    return contig_header;
}

std::vector<std::string> ReferenceGenome::getChromosomes() const
{
    return this->chromosomes;
}

uint32_t ReferenceGenome::getChromosomeLength(std::string chr) const
{
    try
    {
        return this->chr_to_length.at(chr);
    }
    catch (const std::out_of_range& e)
    {
        printError("Chromosome " + chr + " not found in reference genome");
        return 0;
    }
}
