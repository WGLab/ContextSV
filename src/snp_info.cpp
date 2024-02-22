#include "snp_info.h"

/// @cond
#include <string>
#include <tuple>
#include <mutex>
#include <iostream>
/// @endcond

#define MIN_PFB 0.01


// Function to remove the 'chr' prefix from chromosome names
std::string removeChrPrefix(std::string chr)
{
    if (chr.find("chr") != std::string::npos) {
        return chr.substr(3);
    }
    return chr;
}

void SNPInfo::insertSNPAlleleFrequency(std::string chr, int64_t pos, double baf)
{
    chr = removeChrPrefix(chr);

    // Add the chromosome to the SNP B-allele frequency map if it does not exist
    if (this->snp_baf_map.find(chr) == this->snp_baf_map.end()) {
        this->snp_baf_map[chr] = BST();
    }

    // Insert the SNP into the map with its position and B-allele frequency
    // using a binary search tree to keep the SNP positions sorted
    this->snp_baf_map[chr].insert({pos, baf});
}

void SNPInfo::insertSNPPopulationFrequency(std::string chr, int64_t pos, double pfb)
{
    chr = removeChrPrefix(chr);

    // Add the chromosome to the SNP population frequency map if it does not
    // exist
    if (this->snp_pfb_map.find(chr) == this->snp_pfb_map.end()) {
        this->snp_pfb_map[chr] = std::unordered_map<int64_t, double>();
    }

    // Insert the SNP into the map with its position and population frequency of
    // the B allele
    this->snp_pfb_map[chr][pos] = pfb;
}

std::tuple<std::vector<int64_t>, std::vector<double>, std::vector<double>> SNPInfo::querySNPs(std::string chr, int64_t start, int64_t end)
{
    // Lock the mutex for reading SNP information
    // std::lock_guard<std::mutex> lock(this->snp_info_mtx);

    chr = removeChrPrefix(chr);

    // Create an ordered map of SNP positions to BAF and PFB values
    std::map<int64_t, std::tuple<double, double>> snp_map;

    // Query SNPs within a range (start, end) and return their BAF and PFB
    // values as separate vectors
    std::vector<double> bafs;
    std::vector<double> pfbs;
    std::vector<int64_t> pos;
    
    // Check if the chromosome exists in the B-allele frequency map
    if (this->snp_baf_map.find(chr) == this->snp_baf_map.end()) {
        std::cerr << "Chromosome " << chr << " not found in SNP BAF map" << std::endl;
        return std::make_tuple(pos, bafs, pfbs);
    }

    // Query the SNPs within the range and return their BAFs and corresponding
    // positions
    auto& baf_bst = this->snp_baf_map[chr];
    auto baf_start = baf_bst.lower_bound({start, 0.0});
    auto baf_end = baf_bst.upper_bound({end, 0.0});
    for (auto it = baf_start; it != baf_end; it++) {
        bafs.push_back(std::get<1>(*it));
        pos.push_back(std::get<0>(*it));
        //std::cout << "SNP at " << std::get<0>(*it) << " with BAF " << std::get<1>(*it) << std::endl;
    }

    // Define a default PFB value for SNPs with no population frequency data
    pfbs = std::vector<double>(bafs.size(), MIN_PFB);

    // Check if the chromosome exists in the population frequency map
    if (this->snp_pfb_map.find(chr) == this->snp_pfb_map.end()) {
        std::cerr << "Chromosome " << chr << " not found in SNP PFB map" << std::endl;
        return std::make_tuple(pos, bafs, pfbs);
    }

    // Query the PFBs for all SNP positions with PFB data
    auto& pfb_map = this->snp_pfb_map[chr];
    for (size_t i = 0; i < pos.size(); i++) {
        int64_t snp_pos = pos[i];
        if (pfb_map.find(snp_pos) != pfb_map.end()) {
            pfbs[i] = pfb_map[snp_pos];
            //std::cout << "SNP at " << snp_pos << " with PFB " << pfb_map[snp_pos] << std::endl;
        }
    }
    
    return std::make_tuple(pos, bafs, pfbs);
}

std::vector<std::string> SNPInfo::getChromosomes()
{
    // Get a list of chromosomes with SNP information
    std::vector<std::string> chromosomes;
    for (auto it = this->snp_baf_map.begin(); it != this->snp_baf_map.end(); it++) {
        chromosomes.push_back(it->first);
    }
    return chromosomes;
}

std::pair<int64_t, int64_t> SNPInfo::getSNPRange(std::string chr)
{
    chr = removeChrPrefix(chr);

    // Get the range of SNP positions for a given chromosome
    int64_t start = 0;
    int64_t end = 0;
    if (this->snp_baf_map.find(chr) != this->snp_baf_map.end()) {
        auto& baf_bst = this->snp_baf_map[chr];
        start = std::get<0>(*baf_bst.begin());
        end = std::get<0>(*baf_bst.rbegin());
    }
    return std::make_pair(start, end);
}
