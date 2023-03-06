

#include "lrr.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define BUFFER_SIZE 1024


LogRRatio::LogRRatio()
= default;

std::vector<int> LogRRatio::getNthReadCoverage(int read_index, std::string input_filepath)
{
    char cmd[BUFFER_SIZE];
    FILE *fp;
    char line[BUFFER_SIZE];
    char chr [20];
    int pos, depth;

    fprintf(stdout, "file = %s\n", input_filepath.c_str());

    // Call samtools depth to calculate coverage
    snprintf(cmd, BUFFER_SIZE, "samtools depth %s", input_filepath.c_str());    
    fp = popen(cmd, "r");
    if (fp == NULL) {
        fprintf(stderr, "Failed to run samtools command\n");
        exit(EXIT_FAILURE);
    }

    // Parse output of samtools depth
    while (fgets(line, BUFFER_SIZE, fp) != NULL) {
        if (sscanf(line, "%s%d%d", chr, &pos, &depth) == 3) {
            printf("%s\t%d\t%d\n", chr, pos, depth);
        }
    }

    // Close file and clean up
    pclose(fp);

    return std::vector<int>();
}
