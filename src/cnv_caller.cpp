

#include "cnv_caller.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>       /* log2 */
#include <htslib/sam.h>

#define BUFFER_SIZE 1024

CNVCaller::CNVCaller() = default;

std::vector<double> CNVCaller::run(std::string input_filepath, SNVCaller snv_obj)
{
    // Get alignment endpoints
    getAlignmentEndpoints(input_filepath);

    // Calculate LRRs
    std::vector<double> log_r_ratios;
    log_r_ratios = calculateLogRRatios(input_filepath);

    return log_r_ratios;
}

std::vector<double> CNVCaller::calculateLogRRatios(std::string input_filepath)
{
    char target_chr[] = "chr1";
    std::vector<double> chr_lrr;
    
    // Calculate mean chromosome coverage
    RegionCoverage chr_cov = getMeanCoverage(input_filepath, target_chr);

    // Calculate Log R ratios
    int start_pos = this->align_start, end_pos = (start_pos + window_size - 1);
    while (end_pos < this->align_end) {

        // Calculate window mean coverage
        RegionCoverage region_cov = getMeanCoverage(input_filepath, target_chr, start_pos, end_pos);
        
        if (region_cov.valid) {
            // Calculate the region's Log R ratio
            double region_lrr = log2(region_cov.mean / chr_cov.mean);
            chr_lrr.push_back(region_lrr);
            fprintf(stdout, "%.3f/%.3f\n", region_cov.mean, chr_cov.mean);
            fprintf(stdout, "LRR=%.3f\n", region_lrr);
        }

        // Update indexes
        start_pos = (end_pos + 1);
        end_pos += window_size;
    }

    std::cout << "CNVs calculated" << std::endl;

    return chr_lrr;
}

RegionCoverage CNVCaller::getMeanCoverage(std::string input_filepath, char* chr, int start_pos, int end_pos)
{
    char cmd[BUFFER_SIZE];
    FILE *fp;
    char line[BUFFER_SIZE];
    RegionCoverage cov;
    bool log_debug = false;  // Log debugging output

    // Open a SAMtools process to calculate cumulative read depth and position
    // counts (non-zero depths only) for a single chromosome
    if (end_pos == 1) {
        // Run the entire chromosome
        log_debug = true;
        snprintf(cmd, BUFFER_SIZE,\
        "samtools depth -r %s %s | awk '{c++;s+=$3}END{print c, s}'",\
        chr, input_filepath.c_str());  // Remove '-a' for debugging
    } else {
        // Run on a region of the chromosome
        snprintf(cmd, BUFFER_SIZE,\
        "samtools depth -r %s:%d-%d %s | awk '{c++;s+=$3}END{print c, s}'",\
        chr, start_pos, end_pos, input_filepath.c_str());  // Remove '-a' for debugging
    }
    
    if (log_debug) {fprintf(stdout, "%s\n", cmd);};  // Print the command
    
    fp = popen(cmd, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Failed to run command\n");
        exit(EXIT_FAILURE);
    }

    // Parse the outputs
    int pos_count, cum_depth;
    if (fgets(line, BUFFER_SIZE, fp) != NULL)
    {
        if (sscanf(line, "%d%d", &pos_count, &cum_depth) == 2)
        {
            // Calculate the mean chromosome coverage
            double chr_mean_coverage = double(cum_depth) / pos_count;
            fprintf(stdout, "%s mean coverage = %.3f\ncumulative depth = %d\nregion length=%d\n", chr, chr_mean_coverage, cum_depth, pos_count);
            cov.length = pos_count;
            cov.mean   = chr_mean_coverage;
            cov.valid  = true;
        }
    }
    pclose(fp);  // Close the process

    return cov;
}

int CNVCaller::getAlignmentEndpoints(std::string input_filepath)
{
	samFile *fp_in = hts_open(input_filepath.c_str(), "r");  // Open BAM file
	bam_hdr_t *bamHdr = sam_hdr_read(fp_in);  // Read header
    bam1_t *aln = bam_init1();  // Initialize alignment record
    std::vector<int> end_positions (2);
	
    int read_count = 0;
    int primary_count = 0;
    int first_position = 0;
	while(sam_read1(fp_in,bamHdr,aln) > 0){
		int32_t pos = aln->core.pos +1; //left most position of alignment in zero-based coordinates (+1)
		char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
		uint32_t len = aln->core.l_qseq; //length of the read
        int map_flag = aln->core.flag;  // Alignment type flag
        if (int(len) > 0) {
            read_count++;

            // Get primary alignments only
            if (! ((map_flag & BAM_FSECONDARY) || (map_flag & BAM_FSUPPLEMENTARY)) )
            {
                primary_count++;
                // Get the first alignment position
                if (std::string(chr) == "chr1" && first_position == 0) {
                        first_position = aln->core.pos + 1;
                        this->align_start = first_position;
                }
            }
        }
	}

    // Get the last alignment position
    int last_position = aln->core.pos + bam_cigar2rlen(aln->core.n_cigar, bam_get_cigar(aln)) - 1;
    this->align_end = last_position;
    std::cout << "alignment start = " << this->align_start << std::endl;
    std::cout << "alignment end   = " << this->align_end << std  ::endl;
	
	bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
	sam_close(fp_in);

    return 0;
}
