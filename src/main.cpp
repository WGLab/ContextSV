
#include "swig_interface.h"

/// @cond DOXYGEN_IGNORE
#include <iostream>
#include <string>

// For signal handling
#include <signal.h>
#include <execinfo.h>

// #include <optional>
/// @endcond

#include "input_data.h"
#include "version.h"
#include "utils.h"


void printStackTrace(int sig)
{
    void *array[10];
    size_t size;

    // get void*'s for all entries on the stack
    size = backtrace(array, 10);

    // print out all the frames to stderr
    fprintf(stderr, "Error: signal %d:\n", sig);
    backtrace_symbols_fd(array, size, STDERR_FILENO);
    exit(1);
}


void runContextSV(const std::unordered_map<std::string, std::string>& args)
{
    // Set up signal handling
    signal(SIGSEGV, printStackTrace);
    signal(SIGABRT, printStackTrace);
    signal(SIGINT, printStackTrace);
    signal(SIGTERM, printStackTrace);
    signal(SIGILL, printStackTrace);
    signal(SIGFPE, printStackTrace);
    signal(SIGBUS, printStackTrace);

    std::cout << R"(`
    ___         _           _   _____   __
    / __|___ _ _| |_ _____ _| |_/ __\ \ / /
   | (__/ _ \ ' \  _/ -_) \ /  _\__ \\ V / 
    \___\___/_||_\__\___/_\_\\__|___/ \_/  
                                             
    )" << std::endl;

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

    // DBSCAN parameters
    if (args.find("epsilon") != args.end()) {
        input_data.setDBSCAN_Epsilon(std::stod(args.at("epsilon")));
    }

    if (args.find("min-pts") != args.end()) {
        input_data.setDBSCAN_MinPts(std::stoi(args.at("min-pts")));
    }

    if (args.find("min-pts-pct") != args.end()) {
        input_data.setDBSCAN_MinPtsPct(std::stod(args.at("min-pts-pct")));
    }

    // Set up the CNV JSON file if enabled
    if (input_data.getSaveCNVData()) {
        const std::string output_dir = input_data.getOutputDir();
        std::string json_filepath = output_dir + "/CNVCalls.json";
        int json_file_count = 1;
        while (fileExists(json_filepath)) {
            json_filepath = output_dir + "/CNVCalls_" + std::to_string(json_file_count) + ".json";
            json_file_count++;
        }
        input_data.setCNVOutputFile(json_filepath);
        std::cout << "Saving CNV data to: " << json_filepath << std::endl;
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
                << "     --eps <epsilon>             DBSCAN epsilon\n"
                << "     --min-pts <min_pts>         DBSCAN minimum points\n"
                << "     --min-pts-pct <min_pts_pct> Percentage of mean chr. coverage to use for DBSCAN minimum points\n"
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
        } else if (arg == "--min-reads" && i + 1 < argc) {
            args["min-reads"] = argv[++i];
        } else if (arg == "--eps" && i + 1 < argc) {
            args["epsilon"] = argv[++i];
        } else if (arg == "--min-pts" && i + 1 < argc) {
            args["min-pts"] = argv[++i];
        } else if (arg == "--min-pts-pct" && i + 1 < argc) {
            args["min-pts-pct"] = argv[++i];
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
