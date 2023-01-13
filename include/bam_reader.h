//
// Created by jperdomo on 1/8/2023.
//

#ifndef CONTEXTSV_BAM_READER_H
#define CONTEXTSV_BAM_READER_H


#include <string>
class bam_reader {
	private:

	public:
		bam_reader();
		int read(std::string filepath);
};


#endif //CONTEXTSV_BAM_READER_H
