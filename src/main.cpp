
#include "swig_interface.h"
#include "input_data.h"

/// @cond DOXYGEN_IGNORE
#include <iostream>
#include <string>
/// @endcond

// Placeholder for ContextSV library includes
// #include "ContextSV.h"

void runContextSV(const std::string& bamFile, const std::string& refFile, const std::string& vcfFile, const std::string& outputDir, int threadCount = 1, const std::string& hmmFile = "", int windowSize = 2500, int minCNV = 2500, const std::string& eth = "", const std::string& pfbFile = "")
{
    // Placeholder for setting up input data and running ContextSV
    std::cout << "Running ContextSV with the following files:" << std::endl;
    std::cout << "BAM file: " << bamFile << std::endl;
    std::cout << "Reference file: " << refFile << std::endl;
    std::cout << "VCF file: " << vcfFile << std::endl;
    std::cout << "Thread count: " << threadCount << std::endl;
    std::cout << "Output directory: " << outputDir << std::endl;

    // Set up input data
    InputData input_data;
    input_data.setShortReadBam(bamFile);
    input_data.setLongReadBam(bamFile);
    input_data.setRefGenome(refFile);
    input_data.setSNPFilepath(vcfFile);
    //input_data.setChromosome("21");
    //input_data.setRegion("14486099-14515105");
    input_data.setThreadCount(threadCount);
    input_data.setAlleleFreqFilepaths(pfbFile);
    input_data.setHMMFilepath(hmmFile);
    input_data.setOutputDir(outputDir);
    input_data.saveCNVData(false);
    input_data.setThreadCount(threadCount);
    input_data.setWindowSize(windowSize);
    input_data.setMinCNVLength(minCNV);

    // Run ContextSV
    run(input_data);
}

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <bam_file> <ref_file> <vcf_file> <output_dir> <thread_count>" << std::endl;
        return 1;
    }

    std::string bamFile = argv[1];
    std::string refFile = argv[2];
    std::string vcfFile = argv[3];
    std::string outputDir = argv[4];
    int threadCount = std::stoi(argv[5]);

    std::string hmmFile = "";
    int windowSize = 2500;
    int minCNV = 2500;
    std::string eth = "";
    std::string pfbFile = "";
    if (argc == 11) {
        hmmFile = argv[6];
        windowSize = std::stoi(argv[7]);
        minCNV = std::stoi(argv[8]);
        eth = argv[9];
        pfbFile = argv[10];
    }
    
    runContextSV(bamFile, refFile, vcfFile, outputDir, threadCount, hmmFile, windowSize, minCNV, eth, pfbFile);
    
    //std::string hmmFile = argv[6];
    //int windowSize = std::stoi(argv[7]);
    //int minCNV = std::stoi(argv[8]);
    //std::string eth = argv[9];
    //std::string pfbFile = argv[10];
    
    //runContextSV(bamFile, refFile, vcfFile, outputDir, threadCount, "", 2500, 2500, "", "");

    return 0;
}
