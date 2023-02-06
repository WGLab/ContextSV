//
// Created by jperdomo on 1/8/2023.
//

#include "bam_reader.h"

#include <htslib/sam.h>
#include <iostream>
# include <stdio.h>
#include <string>


bam_reader::bam_reader()
= default;

int bam_reader::read(std::string filepath)
{
	samFile *fp_in = hts_open(filepath.c_str(), "r");  // Open BAM file
	bam_hdr_t *bamHdr = sam_hdr_read(fp_in);  // Read header
    bam1_t *aln = bam_init1();  // Initialize alignment record
	
    int read_count = 0;
    int primary_count = 0;
	while(sam_read1(fp_in,bamHdr,aln) > 0){
		
		int32_t pos = aln->core.pos +1; //left most position of alignment in zero based coordianate (+1)
        
		char *chr = bamHdr->target_name[aln->core.tid] ; //contig name (chromosome)
		uint32_t len = aln->core.l_qseq; //length of the read
        
        int map_flag = aln->core.flag;  // Alignment type flag
		
		uint8_t *q = bam_get_seq(aln); //quality string
		uint32_t q2 = aln->core.qual ; //mapping quality
		
		
		char *qseq = (char *)malloc(len);

		for(int i=0; i< len ; i++){
			qseq[i] = seq_nt16_str[bam_seqi(q,i)]; //gets nucleotide id and converts them into IUPAC id.
		}

        //printf("%s\t%d\t%d\t%s\t%s\t%d\n",chr,pos,len,qseq,q,q2);
        // Note: Read count = 106, thus it includes primary + secondary, supplementary alignments as separate records (see: LRS issue)
        int read_length(len);
        printf("Read length = %d\n", read_length);
        if (read_length > 0) {
            read_count++;

            // Get primary alignments only
            if (! ((map_flag & BAM_FSECONDARY) || (map_flag & BAM_FSUPPLEMENTARY)) )
            {
                primary_count++;
            }
        }
	}
    std::cout << "Read count = " << read_count << std::endl;
    std::cout << "Primary count = " << primary_count << std::endl;
	
	bam_destroy1(aln);
	sam_close(fp_in);

	std::cout << "SAMTools succeeded" << std::endl;

    return 0;
}
