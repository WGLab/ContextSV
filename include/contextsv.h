//
// contextsv.h:
// Main class for ContextSV
//

#ifndef CONTEXTSV_H
#define CONTEXTSV_H

#include "input_data.h"
#include "cnv_data.h"
#include "sv_data.h"


class ContextSV {
	private:
		InputData* input_data;

	public:
		ContextSV(InputData& input_data);

		// Entry point
		int run();
};

#endif  // CONTEXTSV_H
