//
// contextsv.h:
// Main class for ContextSV
//

#ifndef CONTEXTSV_H
#define CONTEXTSV_H

#include "common.h"
#include "cnv_map.h"
#include "sv_map.h"

#include <string>

class ContextSV {
	private:
		Common common;		

	public:
		ContextSV(Common common);

		// Entry point
		int run();

		// Integrate CNV and SV calls
		SVMap integrateCNVs(CNVMap cnv_calls, SVMap sv_calls);
};

#endif  // CONTEXTSV_H
