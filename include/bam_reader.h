//
// bam_reader.h:
// Read and perform operations on BAM files.
//

#ifndef CONTEXTSV_BAM_READER_H
#define CONTEXTSV_BAM_READER_H


#include <string>
class BamReader {
	private:

	public:
		BamReader();

		/// Read the BAM file
		int readFile(std::string filepath);
};


#endif //CONTEXTSV_BAM_READER_H
