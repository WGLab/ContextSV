#include "sv_data.h"

/// @cond
#include <iostream>
#include <fstream>
/// @endcond

void SVData::add(std::string chr, int64_t start, int64_t end, int sv_type, std::string alt_allele, std::string data_type)
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
        // Determine the SV length
        int sv_length = end - start;
        if (sv_type == 0) {
            // Deletion
            sv_length *= -1;
        }

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

    // Write the repeat type INFO line
    output_stream << "##INFO=<ID=REPTYPE,Number=1,Type=String,Description=\"Repeat type of the structural variant\">" << std::endl;
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
        int sv_length = info.sv_length;
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

        // If the SV type is unknown, skip it
        if (sv_type == -1) {
            continue;
        }
    
        // Process by SV type
        std::string ref_allele = ".";
        std::string alt_allele = ".";
        std::string repeat_type = "NA";

        // Deletion
        if (sv_type == 0) {
            // Get the reference allele from the reference genome as well as the
            // previous base preceding the SV
            int preceding_pos = std::max(1, pos-1);  // Make sure the position is not negative
            ref_allele = ref_genome.query(chr, preceding_pos, end);

            // Use the previous base as the alternate allele
            alt_allele = ref_allele.substr(0, 1);

            repeat_type = "CONTRAC";  // Deletion
        
        // Duplications and insertions
        } else if (sv_type == 1 || sv_type == 3) {
            // Use the preceding base as the reference allele
            int preceding_pos = std::max(1, pos-1);  // Make sure the position is not negative
            ref_allele = ref_genome.query(chr, preceding_pos, preceding_pos+1);

            // Use the insertion sequence as the alternate allele
            alt_allele = std::get<3>(candidate);

            // Insert the reference base into the alternate allele
            alt_allele.insert(0, ref_allele);

            // Update the end position to the start position to change from
            // query to reference coordinates
            end = pos;

            if (sv_type == 1) {
                repeat_type = "DUP";  // Duplication
            }
        }
        
        // Get the SV type string
        std::string sv_type_str = this->sv_type_map[sv_type];

        // For now, set the QUAL and FILTER as unknown

        // Set the genotype as unspecified for now (Haven't distinguished b/w homozygous, heterozygous)
        std::string genotype = "./.";

        // Create the INFO string
        std::string info_str = "END=" + std::to_string(end) + ";SVTYPE=" + sv_type_str + ";SVLEN=" + std::to_string(sv_length) + ";DP=" + std::to_string(depth) + ";SVMETHOD=" + sv_method + ";ALN=" + data_type_str + ";REPTYPE=" + repeat_type;

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
