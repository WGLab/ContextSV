#ifndef SNP_INFO_H
#define SNP_INFO_H

#include <unordered_map>
#include <map>
#include <string>
#include <vector>
#include <set>
#include <mutex>

// Define the comparator for the binary search tree by SNP position (first
// element of tuple)
struct SNPCompare {
    bool operator()(const std::tuple<uint32_t, double>& a, const std::tuple<uint32_t, double>& b) const {
        return std::get<0>(a) < std::get<0>(b);
    }
};

// Define the data structure for SNP frequencies sorted by position
using BST = std::set<std::tuple<uint32_t, double>, SNPCompare>;

class SNPInfo {
public:
    SNPInfo() {}

    // Insert a SNP into the map with its position and B-allele frequency
    void insertSNPAlleleFrequency(const std::string& chr, uint32_t pos, double baf);

    // Insert a SNP into the map with its position and population frequency of
    // the B allele
    void insertSNPPopulationFrequency(const std::string& chr, uint32_t pos, double pfb);
    
    // Query SNPs within a range (start, end) and return their BAF and PFB values
    std::tuple<std::vector<uint32_t>, std::vector<double>, std::vector<double>> querySNPs(std::string chr, uint32_t start, uint32_t end);

    // Get the range of SNP positions for a given chromosome
    std::pair<uint32_t, uint32_t> getSNPRange(std::string chr);


private:
    // Mutex for reading SNP information
    std::mutex snp_info_mtx;

    // Define the map of chromosome to SNP B-allele frequency
    std::unordered_map<std::string, BST> snp_baf_map;

    // Define the map of chromosome to SNP population frequency
    std::unordered_map<std::string, std::unordered_map<uint32_t, double>> snp_pfb_map;
};

#endif // SNP_INFO_H
