#include "sv_data.h"
#include "vcf_writer.h"

/// @cond
#include <iostream>
#include <fstream>
/// @endcond

int SVData::add(std::string chr, int64_t start, int64_t end, int sv_type, std::string alt_allele, std::string data_type)
{
    // Check if the SV candidate already exists in the map
    SVCandidate candidate(start, end, alt_allele);
    if (this->sv_calls[chr].find(candidate) != this->sv_calls[chr].end()) {
        // Update the alignment-based support count (+1)
        SVInfo& sv_info = this->sv_calls[chr][candidate];
        sv_info.read_support += 1;

        // Update the SV type if it is unknown
        if (sv_info.sv_type == UNKNOWN) {
            sv_info.sv_type = sv_type;
        }

        // Add the alignment type used to call the SV
        sv_info.data_type.insert(data_type);

        return 0;  // SV call already exists

    // Otherwise, add the SV candidate to the map
    } else {
        // For insertions and duplications, the SV length is the length of the
        // inserted sequence, not including the insertion position
        int sv_length = 0;
        if (sv_type == INS || sv_type == DUP || sv_type == TANDUP) {
            sv_length = end - start;
        } else {
            // For deletions, the SV length is the length of the deletion
            sv_length = end - start + 1;
        }

        // Create a new SVInfo object (SV type, alignment support, read depth, data type, SV length, genotype)
        SVInfo sv_info(sv_type, 1, 0, data_type, sv_length, "./.", 0.0);

        // Add the SV candidate to the map
        this->sv_calls[chr][candidate] = sv_info;

        return 1;  // SV call added
    }
}

void SVData::concatenate(const SVData &sv_data)
{
    // Iterate over the chromosomes in the other SVData object
    for (auto const& chr_sv_calls : sv_data.sv_calls) {
        std::string chr = chr_sv_calls.first;

        // Iterate over the SV calls in the other SVData object
        for (auto const& sv_call : chr_sv_calls.second) {

            // Add the SV call to the map of candidate locations. Since the region
            // is unique (per chromosome), there is no need to check if the SV
            // candidate already exists in the map.
            SVCandidate candidate = sv_call.first;  // (start, end, alt_allele)
            SVInfo info = sv_call.second;  // (sv_type, read_support, data_type, sv_length)
            this->sv_calls[chr][candidate] = info;
        }
    }
}

void SVData::updateSVType(std::string chr, SVCandidate candidate, int sv_type, std::string data_type)
{
    // Update the SV type if it is unknown
    SVInfo& sv_info = this->sv_calls[chr][candidate];
    if (sv_info.sv_type == UNKNOWN) {
        sv_info.sv_type = sv_type;

        // Update the alignment type used to call the SV
        sv_info.data_type.insert(data_type);
    }
}

void SVData::updateGenotype(std::string chr, SVCandidate candidate, std::string genotype)
{
    // Update the SV genotype
    // Check if the chromosome and candidate exist in the map
    if (this->sv_calls.find(chr) == this->sv_calls.end() || this->sv_calls[chr].find(candidate) == this->sv_calls[chr].end()) {
        std::cerr << "Error: SV candidate " << chr << ":" << std::get<0>(candidate) << "-" << std::get<1>(candidate) << " does not exist in the SV data" << std::endl;
        return;
    }

    SVInfo& sv_info = this->sv_calls[chr][candidate];
    sv_info.genotype = genotype;
}

void SVData::updateClippedBaseSupport(std::string chr, int64_t pos)
{
    // Update clipped base support
    std::pair<std::string, int64_t> key(chr, pos);
    if (this->clipped_base_support.find(key) != this->clipped_base_support.end()) {
        // Update the depth
        this->clipped_base_support[key] += 1;
    } else {
        // Add the depth
        this->clipped_base_support[key] = 1;
    }
}

int SVData::getClippedBaseSupport(std::string chr, int64_t pos, int64_t end)
{
    // Clipped base support is the maximum clipped base support at the start
    // and end positions
    int clipped_base_support = 0;
    std::pair<std::string, int64_t> pos_key(chr, pos);

    if (pos == end) {
        // If the start and end positions are the same, then the clipped base
        // support is the same at both positions
        clipped_base_support = this->clipped_base_support[pos_key];

    } else{

        // Otherwise, get the clipped base support at the start and end
        // positions
        int pos_support = 0;
        int end_support = 0;
        std::pair<std::string, int64_t> end_key(chr, end);
        if (this->clipped_base_support.find(pos_key) != this->clipped_base_support.end()) {
            pos_support = this->clipped_base_support[pos_key];
        }
        if (this->clipped_base_support.find(end_key) != this->clipped_base_support.end()) {
            end_support = this->clipped_base_support[end_key];
        }
        clipped_base_support = std::max(pos_support, end_support);
    }
    
    return clipped_base_support;
}

void SVData::saveToVCF(FASTAQuery& ref_genome, std::string output_dir)
{
    // Create a VCF writer
    std::string output_vcf = output_dir + "/output.vcf";
    VcfWriter vcf_writer(output_vcf);

    // Set the sample name
    std::string sample_name = "SAMPLE";

    // Set the header lines
    std::vector<std::string> header_lines = {
        std::string("##reference=") + ref_genome.getFilepath(),
        ref_genome.getContigHeader(),
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">",
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">",
        "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Method used to call the structural variant\">",
        "##INFO=<ID=ALN,Number=1,Type=String,Description=\"Alignment type used to call the structural variant\">",
        "##INFO=<ID=CLIPSUP,Number=1,Type=Integer,Description=\"Clipped base support at the start and end positions\">",
        "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description=\"Number of reads supporting the variant\">",
        "##INFO=<ID=REPTYPE,Number=1,Type=String,Description=\"Repeat type\">",
        "##INFO=<ID=HMM,Number=1,Type=Float,Description=\"HMM likelihood\">",
        "##FILTER=<ID=PASS,Description=\"All filters passed\">",
        "##FILTER=<ID=LowQual,Description=\"Low quality\">",
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">"
    };

    // Write the header lines
    vcf_writer.writeHeader(header_lines);

    // Iterate over the SV calls
    std::cout << "Saving SV calls to " << output_vcf << "..." << std::endl;
    std::string sv_method = "CONTEXTSVv0.1";
    int num_sv_calls = this->totalCalls();
    int skip_count = 0;

    // Iterate over the chromosomes
    std::set<std::string> chrs = this->getChromosomes();
    for (auto const& chr : chrs) {
        // Skip the chromosome if there are no SV calls
        if (this->sv_calls.find(chr) == this->sv_calls.end()) {
            continue;
        }

        // Print the chromosome and number of SV calls
        std::cout << "Saving SV calls for " << chr << " (" << this->sv_calls[chr].size() << " SV calls)..." << std::endl;

        // Iterate over the SV calls
        for (auto const& sv_call : this->sv_calls[chr]) {

            // Get the SV candidate and SV info
            SVCandidate candidate = sv_call.first;
            SVInfo info = sv_call.second;
            int sv_type = info.sv_type;
            int read_support = info.read_support;
            int read_depth = info.read_depth;
            int sv_length = info.sv_length;
            std::set<std::string> data_type = info.data_type;
            std::string genotype = info.genotype;
            double hmm_likelihood = info.hmm_likelihood;

            // Convert the data type set to a string
            std::string data_type_str = "";
            for (auto const& type : data_type) {
                data_type_str += type + ",";
            }

            // Get the CHROM, POS, END, and ALT
            int64_t pos = std::get<0>(candidate);
            int64_t end = std::get<1>(candidate);

            // If the SV type is unknown, skip it
            if (sv_type == UNKNOWN) {
                //std::cerr << "Skipping SV call at " << chr << ":" << pos << "-" <<
                //end << " because the SV type is unknown" << std::endl;
                skip_count += 1;
                continue;
            }

            // Process by SV type
            std::string ref_allele = ".";
            std::string alt_allele = ".";
            std::string repeat_type = "NA";

            // Deletion
            if (sv_type == DEL) {
                // Get the deleted sequence from the reference genome, also including the preceding base
                int64_t preceding_pos = (int64_t) std::max(1, (int) pos-1);  // Make sure the position is not negative
                ref_allele = ref_genome.query(chr, preceding_pos, end);

                // Use the preceding base as the alternate allele 
                if (ref_allele != "") {
                    alt_allele = ref_allele.at(0);
                } else {
                    alt_allele = "<DEL>";  // Use symbolic allele for imprecise deletions
                    std::cerr << "Warning: Reference allele is empty for deletion at " << chr << ":" << pos << "-" << end << std::endl;
                }

                // Make the SV length negative
                sv_length = -1 * sv_length;

                // Update the position
                pos = preceding_pos;

            // Duplications and insertions
            } else if (sv_type == INS || sv_type == DUP || sv_type == TANDUP) {
                // Use the preceding base as the reference allele
                int64_t preceding_pos = (int64_t) std::max(1, (int) pos-1);  // Make sure the position is not negative
                ref_allele = ref_genome.query(chr, preceding_pos, preceding_pos);

                // Format novel insertions
                if (sv_type == INS) {
                    // Use the insertion sequence as the alternate allele
                    alt_allele = std::get<2>(candidate);

                    // Insert the reference base into the alternate allele
                    alt_allele.insert(0, ref_allele);

                    // Update the position
                    pos = preceding_pos;

                    // Update the end position to the start position to change from
                    // query to reference coordinates for insertions
                    end = pos;
                } else if (sv_type == DUP) {
                    // Use a symbolic allele for duplications
                    alt_allele = "<DUP>";

                    // Set the repeat type as an interspersed duplication
                    repeat_type = "INTERSPERSED";
                } else if (sv_type == TANDUP) {
                    // Use a symbolic allele for tandem duplications
                    alt_allele = "<DUP>";

                    // Set the repeat type
                    repeat_type = "TANDEM";
                }
            }

            // Get the clipped base support
            int clipped_base_support = this->getClippedBaseSupport(chr, pos, end);

            // Get the SV type string
            std::string sv_type_str = this->sv_type_map[sv_type];

            // Create the INFO string (TODO: Read depth is currently set to 1,
            // needs implementation to get the actual read depth from the BAM
            // file)
            std::string info_str = "END=" + std::to_string(end) + ";SVTYPE=" + sv_type_str + \
                ";SVLEN=" + std::to_string(sv_length) + ";SUPPORT=" + std::to_string(read_support) + \
                ";SVMETHOD=" + sv_method + ";ALN=" + data_type_str + ";CLIPSUP=" + std::to_string(clipped_base_support) + \
                ";REPTYPE=" + repeat_type + ";HMM=" + std::to_string(hmm_likelihood);
                
            // std::string info_str = "END=" + std::to_string(end) + ";SVTYPE=" + sv_type_str + ";SVLEN=" + std::to_string(sv_length) + ";SUPPORT=" + std::to_string(read_support) + ";SVMETHOD=" + sv_method + ";ALN=" + data_type_str + ";CLIPSUP=" + std::to_string(clipped_base_support);

            // Create the FORMAT string for the sample (GT:DP)
            std::string format_str = "GT:DP";

            // Create the sample string (genotype and read depth)
            std::string sample_str = genotype + ":" + std::to_string(read_depth);
            std::vector<std::string> samples = {sample_str};

            // Write the SV call to the file (CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, SAMPLES)
            vcf_writer.writeRecord(chr, pos, ".", ref_allele, alt_allele, ".", "PASS", info_str, format_str, samples);
        }
    }

    // Print the number of SV calls skipped
    std::cout << "Skipped " << skip_count << " of " << num_sv_calls << " SV calls because the SV type is unknown" << std::endl;

    // Close the output stream
    vcf_writer.close();
}

std::map<SVCandidate, SVInfo>& SVData::getChromosomeSVs(std::string chr)
{
    return this->sv_calls[chr];
}

std::set<std::string> SVData::getChromosomes()
{
    std::set<std::string> chromosomes;
    for (auto const& sv_call : this->sv_calls) {
        chromosomes.insert(sv_call.first);
    }
    return chromosomes;
}

int SVData::totalCalls()
{
    int sv_calls = 0;
    for (auto const& sv_call : this->sv_calls) {
        sv_calls += sv_call.second.size();
    }

    return sv_calls;
}

int SVData::totalDeletions()
{
    int del_calls = 0;
    for (auto const& sv_call : this->sv_calls) {
        for (auto const& sv_info : sv_call.second) {
            if (sv_info.second.sv_type == DEL) {
                del_calls += 1;
            }
        }
    }

    return del_calls;
}

int SVData::totalDeletions(std::string chr)
{
    int del_calls = 0;
    for (auto const& sv_info : this->sv_calls[chr]) {
        if (sv_info.second.sv_type == DEL) {
            del_calls += 1;
        }
    }

    return del_calls;
}
