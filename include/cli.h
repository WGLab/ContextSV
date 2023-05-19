//
// cli.h:
// Command-line interface entrypoint.
//

#ifndef CONTEXTSV_CLI_H
#define CONTEXTSV_CLI_H

#include <htslib/sam.h>
#include <string>

class CLI {
    private:
        std::string output_dir;
        std::string ref_filepath;
        std::string bam_filepath;

    public:
        CLI();

        /// Parse input arguments
        int parse(int argc, char** argv);

        /// Run the CLI
        int run();

		/// Get the command argument input
		static std::string getCmdOption(char ** begin, char ** end, const std::string & option);

        /// Check if the input argument is provided
		static bool cmdOptionExists(char** begin, char** end, const std::string& option);

        /// Check if the filepath exists
        static bool fileExists(const std::string &name);

        std::string get_output_dir();
        std::string get_bam_filepath();
        std::string get_ref_filepath();

        static void printHelpText();

};


#endif //CONTEXTSV_CLI_H
