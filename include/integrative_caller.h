//
// integrative_caller.h:
// Integrate SV calls
//

#ifndef INTEGRATIVE_CALLER_H
#define INTEGRATIVE_CALLER_H

#include "common.h"
#include "cnv_map.h"
#include "sv_map.h"

#include <string>

class IntegrativeCaller {
	private:
		Common common;

	public:
		IntegrativeCaller(Common common);

		// Entry point
		int run();

		// Integrate CNV and SV calls
		SVMap integrateCNVs(CNVMap cnv_calls, SVMap sv_calls);
};

#endif  // INTEGRATIVE_CALLER_H
