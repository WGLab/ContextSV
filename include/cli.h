//
// cli.h:
// Command-line interface entrypoint.
//

#ifndef CLI_H
#define CLI_H

#include "common.h"

#include <string>

int run(std::string bam, std::string snps, std::string outdir, std::string region);

#endif //CLI_H
