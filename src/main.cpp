
#include "swig_interface.h"
#include "input_data.h"

/// @cond DOXYGEN_IGNORE
#include <iostream>
#include <string>
/// @endcond

// Placeholder for ContextSV library includes
// #include "ContextSV.h"

void runContextSV(const std::string& bamFile, const std::string& refFile, const std::string& vcfFile, const std::string& outputDir, int threadCount = 1) {
    // Placeholder for setting up input data and running ContextSV
    std::cout << "Running ContextSV with the following files:" << std::endl;
    std::cout << "BAM file: " << bamFile << std::endl;
    std::cout << "Reference file: " << refFile << std::endl;
    std::cout << "VCF file: " << vcfFile << std::endl;
    std::cout << "Output directory: " << outputDir << std::endl;

    // Set up input data
    InputData input_data;
    input_data.setShortReadBam(bamFile);
    input_data.setLongReadBam(bamFile);
    input_data.setRefGenome(refFile);
    input_data.setSNPFilepath(vcfFile);
    input_data.setChromosome("21");
    input_data.setRegion("14486099-14515105");
    input_data.setThreadCount(1);
    input_data.setAlleleFreqFilepaths("");
    input_data.setHMMFilepath("");
    input_data.setOutputDir(outputDir);
    input_data.saveCNVData(true);
    input_data.setThreadCount(threadCount);

    // Run ContextSV
    run(input_data);
}

int main(int argc, char* argv[]) {
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <bam_file> <ref_file> <vcf_file> <output_dir> <thread_count>" << std::endl;
        return 1;
    }

    std::string bamFile = argv[1];
    std::string refFile = argv[2];
    std::string vcfFile = argv[3];
    std::string outputDir = argv[4];
    int threadCount = std::stoi(argv[5]);
    
    runContextSV(bamFile, refFile, vcfFile, outputDir, threadCount);

    return 0;
}
