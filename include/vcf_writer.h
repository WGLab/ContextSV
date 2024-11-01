/// @cond
#include <string>
#include <vector>
#include <fstream>
/// @endcond

class VcfWriter {
public:
    // Constructor
    VcfWriter(const std::string& filename);
    ~VcfWriter();
    void writeHeader(const std::vector<std::string>& headerLines);
    void writeRecord(const std::string& chrom, int pos, const std::string& id,
                     const std::string& ref, const std::string& alt,
                     const std::string& qual, const std::string& filter,
                     const std::string& info, const std::string& format,
                     const std::vector<std::string>& samples);

private:
    std::ofstream file_stream;
};
