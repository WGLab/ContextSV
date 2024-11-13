#ifndef VCF_WRITER_H
#define VCF_WRITER_H

/// @cond
#include <string>
#include <vector>
#include <fstream>
/// @endcond

class VcfWriter {
public:
    explicit VcfWriter(std::string filename);
    // VcfWriter(const std::string& filename);
    ~VcfWriter();

    // Delete copy constructor and assignment operator
    VcfWriter(const VcfWriter&) = delete;
    VcfWriter& operator=(const VcfWriter&) = delete;

    void writeHeader(const std::vector<std::string>& headerLines);
    void writeRecord(const std::string& chrom, int pos, const std::string& id,
                     const std::string& ref, const std::string& alt,
                     const std::string& qual, const std::string& filter,
                     const std::string& info, const std::string& format,
                     const std::vector<std::string>& samples);

private:
    std::ofstream file_stream;
};

#endif  // VCF_WRITER_H
