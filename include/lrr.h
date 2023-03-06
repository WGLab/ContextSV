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

        /// Get the Nth read coverage
		std::vector<int> getNthReadCoverage(int read_index, std::string input_filepath);
};


#endif //CONTEXTSV_LRR_H

