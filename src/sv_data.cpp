#include "sv_data.h"
#include "vcf_writer.h"

/// @cond
#include <iostream>
#include <fstream>
/// @endcond

void SVData::add(std::string chr, int64_t start, int64_t end, int sv_type, std::string alt_allele, std::string data_type, std::mutex& mtx)
{
    // Lock the mutex for thread safety
    std::lock_guard<std::mutex> lock(mtx);

    // Add the SV call to the map of candidate locations
    SVCandidate candidate(chr, start, end, alt_allele);
    if (this->sv_calls.find(candidate) != this->sv_calls.end()) {

        // Update the read depth
        SVInfo& sv_info = this->sv_calls[candidate];
        sv_info.read_depth += 1;

        // Update the SV type if it is unknown
        if (sv_info.sv_type == -1) {
            sv_info.sv_type = sv_type;
        }

        // Update the alignment type used to call the SV
        sv_info.data_type.insert(data_type);

    } else {
        // Determine the SV length
        int sv_length = end - start + 1;

        // Create a new SVInfo object
        SVInfo sv_info(sv_type, 1, data_type, sv_length);

        // Add the SV candidate to the map
        this->sv_calls[candidate] = sv_info;
    }
}

SVData::SVData(FASTAQuery &ref_genome)
{
    // Set the reference genome
    this->ref_genome = &ref_genome;
}

std::string SVData::getRefGenome()
{
    return this->ref_genome->getFilepath();
}

std::string SVData::getSequence(std::string chr, int64_t pos_start, int64_t pos_end)
{
    // Query the reference genome
    return this->ref_genome->query(chr, pos_start, pos_end);
}

void SVData::updateSVType(SVCandidate candidate, int sv_type, std::string data_type)
{
    // Update the SV type
    SVInfo& sv_info = this->sv_calls[candidate];
    sv_info.sv_type = sv_type;

    // Update the alignment type used to call the SV
    sv_info.data_type.insert(data_type);
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
    std::string output_vcf = output_dir + "/sv_calls.vcf";
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
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">",
        "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Method used to call the structural variant\">",
        "##INFO=<ID=ALN,Number=1,Type=String,Description=\"Alignment type used to call the structural variant\">",
        "##INFO=<ID=REPTYPE,Number=1,Type=String,Description=\"Repeat type of the structural variant\">",
        "##INFO=<ID=CLIPDP,Number=1,Type=Integer,Description=\"Clipped base support at the start and end positions\">",
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
    int num_sv_calls = this->sv_calls.size();
    int skip_count = 0;
    for (auto const& sv_call : this->sv_calls) {

        // Get the SV candidate and SV info
        SVCandidate candidate = sv_call.first;
        SVInfo info = sv_call.second;
        int sv_type = info.sv_type;
        int depth = info.read_depth;
        int sv_length = info.sv_length;
        std::set<std::string> data_type = info.data_type;

        // Convert the data type set to a string
        std::string data_type_str = "";
        for (auto const& type : data_type) {
            data_type_str += type + ",";
        }

        // Get the CHROM, POS, END, and ALT
        std::string chr = std::get<0>(candidate);
        int64_t pos = std::get<1>(candidate);
        int64_t end = std::get<2>(candidate);

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
            // Get the reference allele from the reference genome as well as the
            // previous base preceding the SV
            int64_t preceding_pos = (int64_t) std::max(1, (int) pos-1);  // Make sure the position is not negative
            ref_allele = ref_genome.query(chr, preceding_pos, end);

            // Use the previous base as the alternate allele
            alt_allele = ref_allele.substr(0, 1);

            // Make the SV length negative
            sv_length = -1 * sv_length;
            
            // Update the position
            pos = preceding_pos;

            repeat_type = "CONTRAC";  // Deletion
        
        // Duplications and insertions
        } else if (sv_type == INS || sv_type == DUP) {
            // Use the preceding base as the reference allele
            int64_t preceding_pos = (int64_t) std::max(1, (int) pos-1);  // Make sure the position is not negative
            ref_allele = ref_genome.query(chr, preceding_pos, preceding_pos);

            // Use the insertion sequence as the alternate allele
            alt_allele = std::get<3>(candidate);

            // Insert the reference base into the alternate allele
            alt_allele.insert(0, ref_allele);
            
            // Update the position
            pos = preceding_pos;

            // Update the end position to the start position to change from
            // query to reference coordinates
            end = pos;

            if (sv_type == DUP) {
                repeat_type = "DUP";  // Duplication
            }
        }

        // Get the clipped base support
        int clipped_base_support = this->getClippedBaseSupport(chr, pos, end);
        
        // Get the SV type string
        std::string sv_type_str = this->sv_type_map[sv_type];

        // Set the genotype as unspecified for now (Haven't distinguished b/w homozygous, heterozygous)
        std::string genotype = "./.";

        // Create the INFO string
        std::string info_str = "END=" + std::to_string(end) + ";SVTYPE=" + sv_type_str + ";SVLEN=" + std::to_string(sv_length) + ";DP=" + std::to_string(depth) + ";SVMETHOD=" + sv_method + ";ALN=" + data_type_str + ";REPTYPE=" + repeat_type + ";CLIPDP=" + std::to_string(clipped_base_support);
        //std::string info_str = "END=" + std::to_string(end) + ";SVTYPE=" + sv_type_str + ";SVLEN=" + std::to_string(sv_length) + ";DP=" + std::to_string(depth) + ";SVMETHOD=" + sv_method + ";ALN=" + data_type_str + ";REPTYPE=" + repeat_type;

        // Create the FORMAT string
        std::string format_str = "GT:DP";

        // Create the sample string
        std::string sample_str = genotype + ":" + std::to_string(depth);
        std::vector<std::string> samples = {sample_str};

        // Write the SV call to the file
        vcf_writer.writeRecord(chr, pos, ".", ref_allele, alt_allele, ".", ".", info_str, format_str, samples);
        //output_stream << chr << "\t" << pos << "\t" << "." << "\t" << ref_allele << "\t" << alt_allele << "\t" << "." << "\t" << "." << "\t" << info_str << "\t" << format_str << "\t" << sample_str << std::endl;
    }

    // Print the number of SV calls skipped
    std::cout << "Skipped " << skip_count << " of " << num_sv_calls << " SV calls because the SV type is unknown" << std::endl;

    // Close the output stream
    vcf_writer.close();
    // output_stream.close();
    std::cout << "Saved SV calls to " << output_vcf << std::endl;
}
