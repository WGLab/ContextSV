// SVCaller: Detect SVs from long read alignments

#ifndef SV_CALLER_H
#define SV_CALLER_H

#include "input_data.h"
#include "cnv_data.h"
#include "sv_data.h"

// Type aliases
class SVCaller;
        // Alignment location (chr, start, end)
        typedef std::tuple<std::string, int, int> AlignmentLocation;

        // Alignment vector (alignment location)
        typedef std::vector<AlignmentLocation> AlignmentVector;

        // Query map (query name, alignment vector)
        typedef std::map<std::string, AlignmentVector> QueryMap;

class SVCaller {
    private:
        //int max_indel_dist = 1000;  // Maximum distance between two indels to
        //be considered as a single SV
        int max_indel_dist = 10;  // Maximum distance between two indels to be considered as a single SV
        //int min_sv_size = 50;       // Minimum SV size to be considered
        //int min_sv_size = 30;       // Minimum SV size to be considered
        int min_sv_size = 50;       // Minimum SV size to be considered
        int min_mapq = 20;          // Minimum mapping quality to be considered
        InputData input_data;

        // Detect SVs from long read alignments in the CIGAR string
        SVData detectSVsFromCIGAR(std::string chr, int32_t pos, uint32_t* cigar, int cigar_len);

        // Detect SVs from split-read alignments (primary and supplementary)
        SVData detectSVsFromSplitReads();

    public:
        SVCaller(InputData input_data);

        // Detect SVs and predict SV type from long read alignments and CNV calls
        SVData run();
};

#endif // SV_CALLER_H
