#ifndef QUERY_MAP_H
#define QUERY_MAP_H

#include <string>
#include <map>
#include <vector>

class QueryMap {
    private:
        // Map of alignments by query name
        // Key: Query name
        // Value: Chromosome, start, end
        std::map<std::string, std::vector<std::tuple<std::string, int, int>>> primary_alignments;
        std::map<std::string, std::vector<std::tuple<std::string, int, int>>> supplementary_alignments;

    public:
        // Add an alignment to the map (0 = primary, 1 = supplementary)
        void addAlignment(std::string chr, int start, int end, int alignment_type);
        std::map<std::string, std::vector<std::tuple<std::string, int, int>>> getPrimaryAlignments();
        std::vector<std::tuple<std::string, int, int>> getSupplementaryAlignments(std::string query_name);
};

#endif // QUERY_MAP_H
