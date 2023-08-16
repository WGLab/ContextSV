#ifndef SV_DATA_H
#define SV_DATA_H

#include <string>
#include <map>


// Type aliases
class SVData;
        // SV candidate location (chr, start, end)
        typedef std::tuple<std::string, int, int> SVCandidate;

        // SV candidate info (SV type, read depth, reference allele, alternate allele)
        typedef std::tuple<int, int, std::string, std::string> SVInfo;

        // SV info map (SV candidate location, SV info)
        typedef std::map<SVCandidate, SVInfo> SVInfoMap;

class SVData {
    private:
        // SV candidate info map
        SVInfoMap sv_calls;

    public:
        void addSVCall(std::string, int start, int end, int sv_type);
        void addSVCalls(SVData sv_calls);

        // Update the SV type for a given SV candidate
        void updateSVType(std::string chr, int start, int end, int sv_type);
        
        // Save SV calls to VCF
        void saveToVCF(std::string output_dir);

        // Begin and end iterators for the SV candidate map
        SVInfoMap::iterator begin() { return this->sv_calls.begin(); }
        SVInfoMap::iterator end() { return this->sv_calls.end(); }
};

#endif // SV_DATA_H
