#include "sv_data.h"

#include <iostream>
#include <fstream>

void SVData::addSVCall(std::string chr, int start, int end, int sv_type, std::string alt_allele)
{
    // TODO: Only need to compare alt. allele in the map's key
    // Add the SV call to the map of candidate locations
    SVCandidate key(chr, start, end, sv_type, alt_allele);
    if (this->sv_calls.find(key) != this->sv_calls.end()) {
        // Update the depth for the SV candidate
        this->sv_calls[key] += 1;
    } else {
        // Add the SV candidate to the map
        this->sv_calls[key] = 1;
    }
}

SVData::SVData(FASTAQuery &ref_genome)
{
    // Set the reference genome
    this->ref_genome = &ref_genome;
}

std::string SVData::getRefGenome()
{
    return this->ref_genome->getFilepath();
}

std::string SVData::getSequence(std::string chr, int pos_start, int pos_end)
{
    // Query the reference genome
    return this->ref_genome->query(chr, pos_start, pos_end);
}

void SVData::updateSVType(SVCandidate key, int sv_type)
{
    // Update the SV type for a given SV candidate
    if (this->sv_calls.find(key) != this->sv_calls.end()) {
        // Update the SV type
        std::get<3>(key) = sv_type;
    } else {
        std::cerr << "Error: SV candidate not found in SVData::updateSVType()" << std::endl;
        exit(1);
    }
}

void SVData::saveToVCF(FASTAQuery& ref_genome, std::string output_dir)
{

    // Create a VCF file for the SV calls
    std::string output_vcf = output_dir + "/sv_calls.vcf";

    // Remove the file if it already exists
    std::remove(output_vcf.c_str());

    // Open the output stream
    std::ofstream output_stream(output_vcf);
    if (!output_stream.is_open()) {
        std::cerr << "Error: Unable to open " << output_vcf << std::endl;
        exit(1);
    }

    // Write the header
    output_stream << "##fileformat=VCFv4.2" << std::endl;
    output_stream << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << std::endl;
    output_stream << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << std::endl;
    output_stream << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">" << std::endl;
    output_stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;

    // Iterate over the SV calls
    std::cout << "Saving SV calls to " << output_vcf << "..." << std::endl;
    for (auto const& sv_call : this->sv_calls) {

        // Get the SV candidate and SV info
        SVCandidate key = sv_call.first;
        int depth = sv_call.second;

        // Get the CHROM, POS, END, SVTYPE, and ALT
        std::string chr = std::get<0>(key);
        int pos = std::get<1>(key);
        int end = std::get<2>(key);
        int sv_type = std::get<3>(key);
        std::string alt_allele = std::get<4>(key);

        // // Get the depth and reference allele
        // int depth = std::get<1>(value);
        // std::string ref_allele = std::get<0>(value);

        // Get the reference allele from the reference genome
        std::cout << "Getting reference allele..." << std::endl;
        //std::string ref_allele = this->getSequence(chr, pos, end);
        std::string ref_allele = ref_genome.query(chr, pos, end);
        std::cout << "ref_allele first 10: " << ref_allele.substr(0, 10) << std::endl;

        // Get the SV length
        int sv_len = end - pos;

        std::cout << "chr: " << chr << ", pos: " << pos << ", end: " << end << ", sv_type: " << sv_type << ", alt_allele: " << alt_allele << ", depth: " << depth << ", ref_allele: " << ref_allele << ", sv_len: " << sv_len << std::endl;

        // Write the SV call to the file
        output_stream << chr << "\t" << pos << "\t" << "." << "\t" << ref_allele << "\t" << alt_allele << "\t" << "." << "\t" << "." << "\t" << "SVTYPE=" << sv_type << ";SVLEN=" << sv_len << ";DP=" << depth << std::endl;
    }

    // Close the output stream
    output_stream.close();
    std::cout << "Saved SV calls to " << output_vcf << std::endl;
}
