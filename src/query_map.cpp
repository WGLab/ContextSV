#include "query_map.h"

void QueryMap::addAlignment(std::string chr, int start, int end, int alignment_type)
{
    std::tuple<std::string, int, int> alignment = std::make_tuple(chr, start, end);
    if (alignment_type == 0) {
        this->primary_alignments["query_name"].push_back(alignment);
    } else {
        this->supplementary_alignments["query_name"].push_back(alignment);
    }
}

std::map<std::string, std::vector<std::tuple<std::string, int, int>>> QueryMap::getPrimaryAlignments()
{
    return this->primary_alignments;
}

std::vector<std::tuple<std::string, int, int>> QueryMap::getSupplementaryAlignments(std::string query_name)
{
    return this->supplementary_alignments[query_name];
}
