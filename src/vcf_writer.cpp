#include "vcf_writer.h"

/// @cond
#include <iostream>
#include <fstream>
/// @endcond

VcfWriter::VcfWriter(const std::string &filename)
{
    // Remove the file if it already exists
    std::remove(filename.c_str());

    // Open the VCF file
    this->file_stream.open(filename);
    if (!this->file_stream.is_open()) {
        std::cerr << "Error: Unable to open " << filename << std::endl;
        exit(1);
    }
}

void VcfWriter::writeHeader(const std::vector<std::string> &headerLines)
{
    // Add the file format
    std::string file_format = "##fileformat=VCFv4.2";
    this->file_stream << file_format << std::endl;

    // Add date and time
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];
    time (&rawtime);
    timeinfo = localtime(&rawtime);
    strftime(buffer, sizeof(buffer), "%Y%m%d", timeinfo);
    file_stream << "##fileDate=" << buffer << std::endl;

    // Add source
    std::string source = "##source=ContexSV";
    this->file_stream << source << std::endl;

    // Loop over the header metadata lines
    for (auto &line : headerLines) {
        this->file_stream << line << std::endl;
    }

    // Add the header line
    std::string header_line = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE";
    this->file_stream << header_line << std::endl;

    // Flush the stream to ensure that the header is written
    this->file_stream.flush();
}

void VcfWriter::writeRecord(const std::string &chrom, int pos, const std::string &id, const std::string &ref, const std::string &alt, const std::string &qual, const std::string &filter, const std::string &info, const std::string &format, const std::vector<std::string> &samples)
{
    // Write a record to the VCF file
    this->file_stream << chrom << "\t" << pos << "\t" << id << "\t" << ref << "\t" << alt << "\t" << qual << "\t" << filter << "\t" << info << "\t" << format << "\t" << samples[0] << std::endl;
}

void VcfWriter::close()
{
    // Close the VCF file
    this->file_stream.close();
}
