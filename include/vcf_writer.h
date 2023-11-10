/// @cond
#include <string>
#include <vector>
#include <fstream>
/// @endcond

class VcfWriter {
public:
    // Constructor
    VcfWriter(const std::string& filename);
    void writeHeader(const std::vector<std::string>& headerLines);
    void writeRecord(const std::string& chrom, int pos, const std::string& id,
                     const std::string& ref, const std::string& alt,
                     const std::string& qual, const std::string& filter,
                     const std::string& info, const std::string& format,
                     const std::vector<std::string>& samples);

    // Close the VCF file
    void close();

private:
    std::ofstream file_stream;
};
