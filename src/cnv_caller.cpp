
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

CNVCaller::CNVCaller(Common common, std::vector<int> snp_positions)
{
    this->common = common;
    this->snp_positions = snp_positions;
}

std::vector<double> CNVCaller::run()
{
    // Get alignment endpoints
    // std::cout << "Getting alignment endpoints..." << std::endl;
    // getAlignmentEndpoints();
    // std::cout << "Alignment endpoints retrieved." << std::endl;

    // Calculate LRRs
    std::cout << "Calculating LRRs..." << std::endl;
    std::vector<double> log_r_ratios;
    log_r_ratios = calculateLogRRatios();
    std::cout << "LRRs calculated." << std::endl;

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

std::vector<double> CNVCaller::calculateLogRRatios()
{
    std::string target_chr = this->common.get_region_chr();
    int window_size = this->common.get_window_size();
    std::vector<double> chr_lrr;
    
    // Calculate mean chromosome coverage
    std::string input_filepath = this->common.get_bam_filepath();
    fprintf(stdout, "\nCalculating coverage for chromosome: %s\n", target_chr.c_str());
    RegionCoverage chr_cov = getChromosomeCoverage();

    // Set up the output LRR CSV
    std::ofstream lrr_output;
    std::string output_dir = this->common.get_output_dir();
    std::string output_filepath = output_dir + "/lrr_output.csv";
    lrr_output.open(output_filepath);
    lrr_output << "start,end,lrr,baf\n";  // Write headers

    // Set up the output BAF CSV
    std::ofstream baf_output;
    std::string baf_output_filepath = output_dir + "/baf_output.csv";
    baf_output.open(baf_output_filepath);
    baf_output << "pos,baf\n";  // Write headers

    // Determine the region to analyze
    int region_start = this->common.get_region_start();
    int region_end = this->common.get_region_end();
    if (region_start == -1 || region_end == -1) {
        // No region was specified, so analyze the whole chromosome
        region_start = 0;
        region_end = chr_cov.length;
        // region_start = this->align_start;
        // region_end = this->align_end;
    }

    // Loop through the regions
    //int start_pos = this->align_start, end_pos = (start_pos + window_size -
    //1);
    int start_pos = region_start, end_pos = std::min(start_pos + window_size - 1, region_end);
    while (end_pos <= region_end) {

        // Calculate window mean coverage
        // fprintf(stdout, "\nCalculating coverage for region %s:%d-%d\n", target_chr, start_pos, end_pos);
        RegionCoverage region_cov = getRegionCoverage(start_pos, end_pos);
        
        if (region_cov.valid) {
            // Calculate the region's Log R ratio
            double region_lrr = log2(region_cov.mean / chr_cov.mean);
            chr_lrr.push_back(region_lrr);
            fprintf(stdout, "%.3f/%.3f\n", region_cov.mean, chr_cov.mean);
            fprintf(stdout, "LRR=%.3f\n", region_lrr);

            // Save region start, end, LRR, and BAF to CSV
            lrr_output << start_pos << "," << end_pos << "," << region_lrr << "," << region_cov.baf << "\n";

            // Save BAF to CSV from the region's BAF vector
            for (int i = 0; i < region_cov.baf_by_pos.size(); i++)
            {
                // Get the position and BAF
                int pos = region_cov.baf_by_pos[i].first;
                double baf = region_cov.baf_by_pos[i].second;

                // Save to CSV
                baf_output << pos << "," << baf << "\n";
            }
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
RegionCoverage CNVCaller::getChromosomeCoverage()
{
    std::string chr = this->common.get_region_chr();
    std::string input_filepath = this->common.get_bam_filepath();

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
    chr.c_str(), input_filepath.c_str());  // Remove '-a' for debugging

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
            fprintf(stdout, "%s mean coverage = %.3f\ncumulative depth = %d\nregion length=%d\n", chr.c_str(), chr_mean_coverage, cum_depth, pos_count);
            cov.length = pos_count;
            cov.mean   = chr_mean_coverage;
            cov.valid  = true;  // Set the region as valid

            // Print values
            fprintf(stdout, "pos_count = %d\ncum_depth = %d\nmean_coverage (%s) = %.3f\n", pos_count, cum_depth, chr.c_str(), chr_mean_coverage);
        }
    }
    pclose(fp);  // Close the process

    return cov;
}

RegionCoverage CNVCaller::getRegionCoverage(int start_pos, int end_pos)
{
    char cmd[BUFFER_SIZE];
    FILE *fp;
    char line[BUFFER_SIZE];
    RegionCoverage cov;
    bool log_debug = false;  // Log debugging output

    // Get the path to the reference genome
    std::string ref_genome_path = this->common.get_ref_filepath();

    // Run on a region of the chromosome using mpileup
    std::string chr = this->common.get_region_chr();
    std::string input_filepath = this->common.get_bam_filepath();
    snprintf(cmd, BUFFER_SIZE,\
    "samtools mpileup -r %s:%d-%d -f %s %s --no-output-ins --no-output-del --no-output-ends",\
    chr.c_str(), start_pos, end_pos, ref_genome_path.c_str(), input_filepath.c_str());
    
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

    // Vector for storing BAFs and corresponding positions
    std::vector<std::pair<int, double>> baf_pos;

    // Parse the mpileup output
    int pos_count = 0;
    int cum_depth = 0;
    int cum_non_ref = 0;
    double mean_baf = 0.0;
    int snp_count = 0;
    int pos = 0;
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
            // Get the position from column 2
            if (col == 1)
            {
                pos = atoi(tok);
            }

            // Get the reference base from column 3
            else if (col == 2)
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
                // (= SNP)
                double baf = 0.0;
                if (non_ref_count > 0)
                {
                    snp_count++;
                    baf = double(non_ref_count) / (ref_count + non_ref_count);
                    //std::cout << "BAF = " << baf << std::endl;


                    // Store the BAF and its position
                    baf_pos.push_back(std::make_pair(pos, baf));


                    // Update the mean BAF
                    mean_baf = (mean_baf * (snp_count - 1) + baf) / snp_count;
                    //mean_baf = (mean_baf * (pos_count - 1) + baf) / pos_count;
                }
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
        cov.length     = pos_count;
        cov.mean       = chr_mean_coverage;
        cov.baf        = mean_baf;
        cov.baf_by_pos = baf_pos;
    }

    return cov;
}

int CNVCaller::getAlignmentEndpoints()
{
    // TODO: This function only works for the first chromosome in the BAM file.
    // Update to work for any chromosome or remove it.

    // Open the BAM file
    std::string input_filepath = this->common.get_bam_filepath();
    std::cout << "Opening bam file..." << std::endl;
	samFile *fp_in = hts_open(input_filepath.c_str(), "r");  // Open BAM file

    // Check if the file was opened
    if (fp_in == NULL) {
        std::cerr << "Error: failed to open " << input_filepath << "\n";
        return 1;
    }

    std::cout << "Opened." << std::endl;
    std::cout << "Reading header..." << std::endl;
	bam_hdr_t *bamHdr = sam_hdr_read(fp_in);  // Read header
    std::cout << "Header read." << std::endl;

    bam1_t *aln = bam_init1();  // Initialize alignment record
    std::vector<int> end_positions (2);
	
    int read_count = 0;
    int primary_count = 0;
    int first_position = 0;
    int invalid_count = 0;  // Number of invalid alignments

    std::cout << "Reading alignments..." << std::endl;

	while(sam_read1(fp_in,bamHdr,aln) > 0){

        // Check if the alignment is valid
        if (aln->core.tid < 0) {
            //std::cerr << "Error: invalid alignment" << std::endl;
            invalid_count++;
            continue;
        }

		int32_t pos = aln->core.pos +1; //left most position of alignment in zero-based coordinates (+1)
		char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
		uint32_t len = aln->core.l_qseq; //length of the read
        int map_flag = aln->core.flag;  // Alignment type flag
        if (int(len) > 0) {
            read_count++;
            //std::cout << "Read count = " << read_count << std::endl;

            // Get primary alignments only
            if (! ((map_flag & BAM_FSECONDARY) || (map_flag & BAM_FSUPPLEMENTARY)) )
            {
                primary_count++;
                //std::cout << "Primary count = " << primary_count << std::endl;

                // Get the first alignment position
                if (first_position == 0) {
                    std::cout << "First position = " << first_position << std::endl;
                    first_position = aln->core.pos + 1;
                    this->align_start = first_position;

                    // End the loop if the first position is found
                }
            }
        }
	}
	
	bam_destroy1(aln);
    bam_hdr_destroy(bamHdr);
	sam_close(fp_in);

    return 0;
}
