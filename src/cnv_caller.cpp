
// #include "khmm.h"
#include "cnv_caller.h"

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>       /* log2 */
#include <htslib/sam.h>

#define BUFFER_SIZE 1024

CNVCaller::CNVCaller() = default;

std::vector<double> CNVCaller::run()
{
    // Get alignment endpoints
    getAlignmentEndpoints(this->bam_filepath);

    // Calculate LRRs
    std::vector<double> log_r_ratios;
    log_r_ratios = calculateLogRRatios(this->bam_filepath);

    // Calculate BAFs
    //std::vector<double> b_allele_freqs;
    //b_allele_freqs = calculateBAFs(input_filepath);

    // Read the HMM from file
    //std::string hmm_filepath = "data/wgs.hmm";
    //CHMM hmm = ReadCHMM(hmm_filepath.c_str());

    // Set up the input variables
    int num_probes = log_r_ratios.size();
    double *lrr = &log_r_ratios[0];
    double *baf = NULL;
    double *pfb = NULL;
    int *snpdist = NULL;
    double logprob = 0.0;

    // Run the Viterbi algorithm
    //testVit_CHMM(hmm, num_probes, lrr, baf, pfb, snpdist, &logprob);

    // Estimate the hidden states from the LRRs
    // TODO: Follow detect_cnv.pl's example and use the Viterbi algorithm
    // https://github.com/WGLab/PennCNV/blob/b6d76b58821deea4f6fe9dc3c241215f25a7fd67/detect_cnv.pl#LL903C20-L903C20

    // #generate CNV calls
	// 			my $probe_count = scalar (@$lrr)-1;
	// 			khmm::testVit_CHMM ($hmm_model, $probe_count, $lrr, $baf, $pfb, $snpdist, \$logprob);
	// 			analyzeStateSequence ($curcnvcall, $curchr, $pfb, $name, $pos, $sample_sex);

    
    //testVit_CHMM(hmm, log_r_ratios);


    return log_r_ratios;
}

std::vector<double> CNVCaller::calculateLogRRatios(std::string input_filepath)
{
    // char target_chr[] = "chr1";
    char target_chr[10];
    if (this->uses_chr_prefix) {
        strcpy(target_chr, "chr1");
    } else {
        strcpy(target_chr, "1");
    }

    std::vector<double> chr_lrr;
    
    // Calculate mean chromosome coverage
    fprintf(stdout, "\nCalculating coverage for chromosome: %s\n", target_chr);
    RegionCoverage chr_cov = getChromosomeCoverage(input_filepath, target_chr);

    // Set up the output CSV
    std::ofstream lrr_output;
    lrr_output.open ("lrr_output.csv");

    // Write headers
    lrr_output << "start,end,lrr,baf\n";

    // Calculate Log R ratios
    int start_pos = this->align_start, end_pos = (start_pos + window_size - 1);
    while (end_pos < this->align_end) {

        // Calculate window mean coverage
        // fprintf(stdout, "\nCalculating coverage for region %s:%d-%d\n", target_chr, start_pos, end_pos);
        RegionCoverage region_cov = getRegionCoverage(input_filepath, target_chr, start_pos, end_pos);
        
        if (region_cov.valid) {
            // Calculate the region's Log R ratio
            double region_lrr = log2(region_cov.mean / chr_cov.mean);
            chr_lrr.push_back(region_lrr);
            fprintf(stdout, "%.3f/%.3f\n", region_cov.mean, chr_cov.mean);
            fprintf(stdout, "LRR=%.3f\n", region_lrr);

            // Save start, end, LRR, and BAF to CSV
            lrr_output << start_pos << "," << end_pos << "," << region_lrr << "," << region_cov.baf << "\n";
        }

        // Update indexes
        start_pos = (end_pos + 1);
        end_pos += window_size;
    }
    
    lrr_output.close();
    std::cout << "CNVs calculated" << std::endl;

    return chr_lrr;
}

/// Calculate the mean chromosome coverage
RegionCoverage CNVCaller::getChromosomeCoverage(std::string input_filepath, char* chr)
{
    char cmd[BUFFER_SIZE];
    FILE *fp;
    char line[BUFFER_SIZE];
    RegionCoverage cov;
    bool log_debug = false;  // Log debugging output

    // Open a SAMtools process to calculate cumulative read depth and position
    // counts (non-zero depths only) for a single chromosome

    // Run the entire chromosome
    log_debug = true;
    snprintf(cmd, BUFFER_SIZE,\
    "samtools depth -r %s %s | awk '{c++;s+=$3}END{print c, s}'",\
    chr, input_filepath.c_str());  // Remove '-a' for debugging

    // Parse the output
    fp = popen(cmd, "r");
    if (fp == NULL) {
        fprintf(stderr, "Failed to run command\n");
        exit(1);
    }
    
    fprintf(stdout, "%s\n", cmd);  // Print the command
    
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
            cov.valid  = true;  // Set the region as valid
        }
    }
    pclose(fp);  // Close the process

    return cov;
}

RegionCoverage CNVCaller::getRegionCoverage(std::string input_filepath, char* chr, int start_pos, int end_pos)
{
    char cmd[BUFFER_SIZE];
    FILE *fp;
    char line[BUFFER_SIZE];
    RegionCoverage cov;
    bool log_debug = false;  // Log debugging output

    // Get the path to the reference genome
    std::string ref_genome_path = this->ref_filepath;

    // Test: "samtools mpileup -r 1 -f data/hs37d5.fa /home/perdomoj/github/SampleData/ContextSV/LargeSV_trimmed.bam --no-output-ins --no-output-del --no-output-ends

    // Run on a region of the chromosome using mpileup
    snprintf(cmd, BUFFER_SIZE,\
    "samtools mpileup -r %s:%d-%d -f %s %s --no-output-ins --no-output-del --no-output-ends",\
    chr, start_pos, end_pos, ref_genome_path.c_str(), input_filepath.c_str());
    
    // Open the process
    fp = popen(cmd, "r");
    if (fp == NULL)
    {
        fprintf(stderr, "Failed to run command\n");
        exit(EXIT_FAILURE);
    }

    // mpileup output columns:
    // 1. chromosome
    // 2. position
    // 3. reference base
    // 4. read depth
    // 5. read bases
    // 6. base qualities
    // 7. mapping qualities
    // 8. read names
    // 9. read flags
    // 10. read tags

    // Parse the mpileup output
    int pos_count = 0;
    int cum_depth = 0;
    int cum_non_ref = 0;
    double mean_baf = 0.0;
    while (fgets(line, BUFFER_SIZE, fp) != NULL)
    {
        // Parse the line
        char *tok = strtok(line, "\t");  // Tokenize the line
        int col = 0;  // Column index
        int depth = 0;  // Read depth at this position
        
        // Skip zero-depth positions
        bool skip_pos = false;

        // Char for the reference base
        char ref_base = 'N';
        while (tok != NULL)
        {

            // Get the reference base from column 3
            if (col == 2)
            {
                ref_base = tok[0];
            }

            // Get the read depth from column 4
            else if (col == 3)
            {
                depth = atoi(tok);
                if (depth > 0)
                {
                    pos_count++;
                    cum_depth += depth;

                    // Reset the skip flag
                    skip_pos = false;
                } else {
                    // Skip this position
                    skip_pos = true;
                }
            }

            // Get the non-reference base count from column 5
            else if ((col == 4) && (!skip_pos))
            {
                // Cast the token to a string
                std::string bases(tok);

                // Counters for reference and non-reference bases
                int ref_count = 0;
                int non_ref_count = 0;

                // Loop through the bases
                for (int i = 0; i < bases.length(); i++)
                {
                    // Get the base
                    char base = bases[i];

                    // Make the base uppercase
                    base = toupper(base);

                    // Check if the base is non-reference by comparing to A, C, G, T
                    if (base == 'A' || base == 'C' || base == 'G' || base == 'T')
                    {
                        non_ref_count++;

                    // Check if the base is the reference base
                    } else if (base == '.') {
                        ref_count++;
                    }
                }

                // Get the BAF if the non-reference count is greater than zero
                double baf = 0.0;
                if (non_ref_count > 0)
                {
                    baf = double(non_ref_count) / (ref_count + non_ref_count);
                    //std::cout << "BAF = " << baf << std::endl;
                }
            
                // Update the mean BAF
                mean_baf = (mean_baf * (pos_count - 1) + baf) / pos_count;
            }

            tok = strtok(NULL, "\t");
            col++;
        }
    }

    // Determine if the region is valid
    if (pos_count == 0)
    {
        cov.valid = false;

    } else {

        cov.valid = true;

        // Calculate the mean chromosome coverage
        double chr_mean_coverage = double(cum_depth) / pos_count;

        // Update the coverage struct
        cov.length = pos_count;
        cov.mean   = chr_mean_coverage;
        cov.baf    = mean_baf;
    }

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
                if (this->uses_chr_prefix == true) {
                    if (std::string(chr) == "chr1" && first_position == 0) {
                        first_position = aln->core.pos + 1;
                        this->align_start = first_position;
                    }
                } else {
                    if (std::string(chr) == "1" && first_position == 0) {
                        first_position = aln->core.pos + 1;
                        this->align_start = first_position;
                    }
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

void CNVCaller::setChrPrefix(bool uses_chr_prefix)
{
    this->uses_chr_prefix = uses_chr_prefix;
}

void CNVCaller::set_bam_filepath(std::string bam_filepath)
{
    this->bam_filepath = bam_filepath;
}

void CNVCaller::set_ref_filepath(std::string ref_filepath)
{
    this->ref_filepath = ref_filepath;
}
