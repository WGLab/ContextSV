//
// Created by jperdomo on 1/8/2023.
//

#ifndef CONTEXTSV_CLI_H
#define CONTEXTSV_CLI_H


#include <string>

class cli {
    private:
        std::string input_filepath;

    public:
        cli();
        void parse(int argc, char** argv);
		static char* getCmdOption(char ** begin, char ** end, const std::string & option);
		static bool cmdOptionExists(char** begin, char** end, const std::string& option);
        std::string getInputFilepath();
        static void printHelpText();
        static bool fileExists(const std::__cxx11::basic_string<char> &name);
};


#endif //CONTEXTSV_CLI_H
