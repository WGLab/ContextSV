//
// cli.h:
// Command-line interface entrypoint.
//

#ifndef CONTEXTSV_CLI_H
#define CONTEXTSV_CLI_H

#include "common.h"

#include <string>

int run(std::string bam, std::string snps, std::string outdir, std::string region);

#endif //CONTEXTSV_CLI_H
