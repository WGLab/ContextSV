//
// contextsv.h:
// Main class for ContextSV
//

#ifndef CONTEXTSV_H
#define CONTEXTSV_H

#include "input_data.h"


class ContextSV {
	public:
		// explicit ContextSV(InputData& input_data);
		ContextSV() = default;

		// Entry point
		int run(const InputData& input_data) const;
};

#endif  // CONTEXTSV_H
