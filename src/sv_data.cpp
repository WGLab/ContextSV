#include "sv_data.h"

/// @cond
#include <iostream>
#include <fstream>
/// @endcond

void SVData::add(std::string chr, int start, int end, int sv_type, std::string alt_allele, std::string data_type)
{
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
        // Create a new SVInfo object
        SVInfo sv_info(sv_type, 1, data_type);

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

std::string SVData::getSequence(std::string chr, int pos_start, int pos_end)
{
    // Query the reference genome
    return this->ref_genome->query(chr, pos_start, pos_end);
}

void SVData::updateSVType(SVCandidate candidate, int sv_type)
{
    // Update the SV type for a given SV candidate
    if (this->sv_calls.find(candidate) != this->sv_calls.end()) {
        // Update the SV type
        SVInfo& sv_info = this->sv_calls[candidate];
        sv_info.sv_type = sv_type;
    } else {
        std::cerr << "Error: Unable to update SV type for SV candidate (" << std::get<0>(candidate) << ", " << std::get<1>(candidate) << ", " << std::get<2>(candidate) << ", " << std::get<3>(candidate) << ")" << std::endl;
    }
}

void SVData::saveToVCF(FASTAQuery& ref_genome, std::string output_dir)
{
    // Create a VCF file for the SV calls
    std::string output_vcf = output_dir + "/sv_calls.vcf";

    // Remove the file if it already exists
    std::cout << "Removing previous VCF..." << std::endl;
    std::remove(output_vcf.c_str());

    // Open the output stream
    std::cout << "Opening VCF..." << std::endl;
    std::ofstream output_stream(output_vcf);
    if (!output_stream.is_open()) {
        std::cerr << "Error: Unable to open " << output_vcf << std::endl;
        exit(1);
    }

    // Set the sample name
    std::string sample_name = "SAMPLE";

    // Write the header
    output_stream << "##fileformat=VCFv4.2" << std::endl;

    // Get the date
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time(&rawtime);
    timeinfo = localtime(&rawtime);

    // Write the date
    strftime(buffer, sizeof(buffer), "%Y%m%d", timeinfo);

    // Write the date and source
    output_stream << "##fileDate=" << buffer << std::endl;
    output_stream << "##source=ContextSV" << std::endl;

    // Write the reference genome
    std::string ref_genome_filepath = ref_genome.getFilepath();
    output_stream << "##reference=" << ref_genome_filepath << std::endl;

    // Write the contig headers
    output_stream << ref_genome.getContigHeader();

    // Write the INFO lines
    output_stream << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">" << std::endl;
    output_stream << "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" << std::endl;
    output_stream << "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">" << std::endl;
    output_stream << "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth at the locus\">" << std::endl;
    output_stream << "##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Method used to call the structural variant\">" << std::endl;

    // Write the alignment type INFO line
    output_stream << "##INFO=<ID=ALN,Number=1,Type=String,Description=\"Alignment type used to call the structural variant\">" << std::endl;
    //output_stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << std::endl;

    // Write the FILTER lines
    output_stream << "##FILTER=<ID=PASS,Description=\"All filters passed\">" << std::endl;
    output_stream << "##FILTER=<ID=LowQual,Description=\"Low quality\">" << std::endl;

    // Write the FORMAT lines
    // Genotype
    output_stream << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;

    // Read depth
    output_stream << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">" << std::endl;

    // Add the header line
    output_stream << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" << sample_name << std::endl;

    // Iterate over the SV calls
    std::cout << "Saving SV calls to " << output_vcf << "..." << std::endl;
    std::string sv_method = "CONTEXTSVv0.1";
    for (auto const& sv_call : this->sv_calls) {

        // Get the SV candidate and SV info
        SVCandidate candidate = sv_call.first;
        SVInfo info = sv_call.second;
        int sv_type = info.sv_type;
        int depth = info.read_depth;
        std::set<std::string> data_type = info.data_type;

        // Convert the data type set to a string
        std::string data_type_str = "";
        for (auto const& type : data_type) {
            data_type_str += type + ",";
        }

        // Get the CHROM, POS, END, and ALT
        std::string chr = std::get<0>(candidate);
        int pos = std::get<1>(candidate);
        int end = std::get<2>(candidate);
        std::string alt_allele = std::get<3>(candidate);

        // // Print if 21:37857255-37857617
        // if (chr == "21" && pos == 37857255 && end == 37857617) {
        //     std::cout << "[saveToVCF] Found 21:37857255-37857617 type: " << sv_type << std::endl;
        // }

        // // Print if 21:42911938-42912387
        // if (chr == "21" && pos == 42911938 && end == 42912387) {
        //     std::cout << "[saveToVCF] Found 21:42911938-42912387 type: " << sv_type << std::endl;
        // }

        // If the SV type is unknown, skip it
        if (sv_type == -1) {
            continue;
        }

        // // Get the depth and reference allele
        // int depth = std::get<1>(value);
        // std::string ref_allele = std::get<0>(value);

        // Get the reference allele from the reference genome as well as the
        // previous base preceding the SV (Positions are already stored as 1-based).
        std::string ref_allele = ref_genome.query(chr, pos-1, end);
    
        // Use the previous base as the alternate allele
        alt_allele = ref_allele.substr(0, 1);

        // Get the SV length
        int sv_len = end - pos;

        // Change the sign of the SV length if it is a deletion
        if (sv_type == 0) {
            sv_len *= -1;
        }

        // Get the SV type
        std::string sv_type_str = this->sv_type_map[sv_type];

        //std::cout << "SV type: " << sv_type_str << std::endl;

        // // Use symbolic ALT alleles for deletions and duplications
        // if (sv_type == 0 || sv_type == 1) {
        //     alt_allele = "<" + sv_type_str + ">";
        // }

        // Set the genotype as unspecified for now (Haven't distinguished b/w homozygous, heterozygous)
        std::string genotype = "./.";

        // Create the INFO string
        std::string info_str = "END=" + std::to_string(end) + ";SVTYPE=" + sv_type_str + ";SVLEN=" + std::to_string(sv_len) + ";DP=" + std::to_string(depth) + ";SVMETHOD=" + sv_method + ";ALN=" + data_type_str;

        // Create the FORMAT string
        std::string format_str = "GT:DP";

        // Create the sample string
        std::string sample_str = genotype + ":" + std::to_string(depth);

        // Write the SV call to the file
        output_stream << chr << "\t" << pos << "\t" << "." << "\t" << ref_allele << "\t" << alt_allele << "\t" << "." << "\t" << "." << "\t" << info_str << "\t" << format_str << "\t" << sample_str << std::endl;
    }

    // Close the output stream
    output_stream.close();
    std::cout << "Saved SV calls to " << output_vcf << std::endl;
}
