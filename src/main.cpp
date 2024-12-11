
#include "swig_interface.h"
#include "input_data.h"
#include "version.h"

/// @cond DOXYGEN_IGNORE
#include <iostream>
#include <string>
// #include <optional>
/// @endcond

// Placeholder for ContextSV library includes
// #include "ContextSV.h"

void runContextSV(const std::unordered_map<std::string, std::string>& args)
{
    // Placeholder for setting up input data and running ContextSV
    std::cout << "ContextSV version " << VERSION << std::endl;
    std::cout << "Input parameters:" << std::endl;
    for (const auto& arg : args) {
        std::cout << arg.first << ": " << arg.second << std::endl;
    }

    // Set up input data
    InputData input_data;
    input_data.setLongReadBam(args.at("bam-file"));
    input_data.setShortReadBam(args.at("bam-file"));
    input_data.setRefGenome(args.at("ref-file"));
    input_data.setSNPFilepath(args.at("snps-file"));
    input_data.setOutputDir(args.at("output-dir"));
    if (args.find("chr") != args.end()) {
        input_data.setChromosome(args.at("chr"));
    }
    if (args.find("region") != args.end()) {
        input_data.setRegion(args.at("region"));
    }
    if (args.find("thread-count") != args.end()) {
        input_data.setThreadCount(std::stoi(args.at("thread-count")));
    }
    if (args.find("hmm-file") != args.end()) {
        input_data.setHMMFilepath(args.at("hmm-file"));
    }
    if (args.find("sample-size") != args.end()) {
        input_data.setSampleSize(std::stoi(args.at("sample-size")));
    }
    if (args.find("min-cnv") != args.end()) {
        input_data.setMinCNVLength(std::stoi(args.at("min-cnv")));
    }
    if (args.find("eth") != args.end()) {
        input_data.setEthnicity(args.at("eth"));
    }
    if (args.find("pfb-file") != args.end()) {
        input_data.setAlleleFreqFilepaths(args.at("pfb-file"));
    }
    if (args.find("save-cnv") != args.end()) {
        input_data.saveCNVData(true);
    }
    if (args.find("debug") != args.end()) {
        input_data.setVerbose(true);
    }

    // Run ContextSV
    run(input_data);
}

void printUsage(const std::string& programName) {
    std::cerr << "Usage: " << programName << " [options]\n"
                << "Options:\n"
                << "  -b, --bam <bam_file>          Long-read BAM file (required)\n"
                << "  -r, --ref <ref_file>          Reference genome FASTA file (required)\n"
                << "  -s, --snp <vcf_file>          SNPs VCF file (required)\n"
                << "  -o, --outdir <output_dir>     Output directory (required)\n"
                << "  -c, --chr <chromosome>        Chromosome\n"
                << "  -r, --region <region>         Region (start-end)\n"
                << "  -t, --threads <thread_count>  Number of threads\n"
                << "  -h, --hmm <hmm_file>          HMM file\n"
                << "  -n, --sample-size <size>      Sample size for HMM predictions\n"
                << "     --min-cnv <min_length>     Minimum CNV length\n"
                << "  -e, --eth <eth_file>          ETH file\n"
                << "  -p, --pfb <pfb_file>          PFB file\n"
                << "     --save-cnv                 Save CNV data\n"
                << "     --debug                    Debug mode with verbose logging\n"
                << "     --version                  Print version and exit\n"
                << "  -h, --help                    Print usage and exit\n";
}

std::unordered_map<std::string, std::string> parseArguments(int argc, char* argv[]) {
    std::unordered_map<std::string, std::string> args;
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];

        // Handle short and long options
        if ((arg == "-b" || arg == "--bam") && i + 1 < argc) {
            args["bam-file"] = argv[++i];
        } else if ((arg == "-r" || arg == "--ref") && i + 1 < argc) {
            args["ref-file"] = argv[++i];
        } else if ((arg == "-s" || arg == "--snp") && i + 1 < argc) {
            args["snps-file"] = argv[++i];
        } else if ((arg == "-o" || arg == "--outdir") && i + 1 < argc) {
            args["output-dir"] = argv[++i];
        } else if ((arg == "-c" || arg == "--chr") && i + 1 < argc) {
            args["chr"] = argv[++i];
        } else if ((arg == "-r" || arg == "--region") && i + 1 < argc) {
            args["region"] = argv[++i];
        } else if ((arg == "-t" || arg == "--threads") && i + 1 < argc) {
            args["thread-count"] = argv[++i];
        } else if ((arg == "-h" || arg == "--hmm") && i + 1 < argc) {
            args["hmm-file"] = argv[++i];
        } else if ((arg == "-n" || arg == "--sample-size") && i + 1 < argc) {
            args["sample-size"] = argv[++i];
        } else if (arg == "--min-cnv" && i + 1 < argc) {
            args["min-cnv"] = argv[++i];
        } else if ((arg == "-e" || arg == "--eth") && i + 1 < argc) {
            args["eth"] = argv[++i];
        } else if ((arg == "-p" || arg == "--pfb") && i + 1 < argc) {
            args["pfb-file"] = argv[++i];
        } else if (arg == "--save-cnv") {
            args["save-cnv"] = "true";
        } else if (arg == "--debug") {
            args["debug"] = "true";
        } else if ((arg == "-v" || arg == "--version")) {
            std::cout << "ContextSV version " << VERSION << std::endl;
            exit(0);
        } else if (arg == "-h" || arg == "--help") {
            printUsage(argv[0]);
            exit(0);
        } else {
            std::cerr << "Unknown option: " << arg << std::endl;
        }
    }

    // Check for required arguments
    bool hasLR = args.find("bam-file") != args.end();
    bool hasOutput = args.find("output-dir") != args.end();
    bool hasRef = args.find("ref-file") != args.end();
    bool hasSNPs = args.find("snps-file") != args.end();
    bool requiredArgs = hasLR && hasOutput && hasRef && hasSNPs;
    if (!requiredArgs) {
        std::cerr << "Missing required argument(s): -b/--bam, -r/--ref, -s/--snp, -o/--outdir" << std::endl;
        exit(1);
    }

    return args;
}

int main(int argc, char* argv[]) {
    auto args = parseArguments(argc, argv);
    runContextSV(args);

    return 0;
}
