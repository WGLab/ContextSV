#include "utils.h"
#include <stdio.h>
#include <string.h>
#include <iostream>


// Print a progress bar
void printProgress(int progress, int total)
{
    // Get the percentage
    float percent = (float)progress / (float)total * 100.0;

    // Get the number of hashes
    int num_hashes = (int)(percent / 2.0);

    // Print the progress bar
    printf("\r[");
    for (int i = 0; i < num_hashes; i++)
    {
        printf("#");
    }
    for (int i = 0; i < 50 - num_hashes; i++)
    {
        printf(" ");
    }
    printf("] %3.2f%%", percent);
    fflush(stdout);

    // Print a new line if finished
    if (progress == total)
    {
        printf("\n");
    }
}