//
// snv_caller.h:
// SNV caller
//

#ifndef CONTEXTSV_SNV_CALLER_H
#define CONTEXTSV_SNV_CALLER_H

#include "common.h"

#include <string>

class SNVCaller {
	private:
		Common common;
		std::vector<int> snp_positions;
	public:
		SNVCaller(Common common);

		int run();
		std::vector<int> get_snp_positions();
};

#endif
