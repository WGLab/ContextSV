//
// LRR.h:
// Functions for calculating the Log R Ratio
//

#ifndef CONTEXTSV_LRR_H
#define CONTEXTSV_LRR_H


#include <string>
#include <vector>

class LogRRatio {
    public:
        LogRRatio();

        /// Calculate the Log R Ratio
		std::vector<int> run(std::string input_filepath);

        /// Calculate per-position coverage using 'samtools depth'
		std::vector<int> getCoverage(std::string input_filepath);
};


#endif //CONTEXTSV_LRR_H

