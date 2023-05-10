//
// integrative_caller.h:
// Integrate SV calls
//

#ifndef CONTEXTSV_INTEGRATIVE_CALLER_H
#define CONTEXTSV_INTEGRATIVE_CALLER_H


#include <string>
class IntegrativeCaller {
	private:
		// Nothing here yet
	public:
		IntegrativeCaller();

		/// Entry point
		int run(std::string filepath);

		/// Check if the bam file uses chr prefix notation
		int bamHasChrPrefix(std::string filepath, bool& uses_chr_prefix);
};

#endif
