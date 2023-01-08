//
// Created by jperdomo on 1/8/2023.
//

#ifndef CONTEXTSV_CLI_H
#define CONTEXTSV_CLI_H


#include <string>

class CLI {
    private:
        std::string input_filepath{};

    public:
        CLI();
        void parse(int argc, char** argv);
        std::string getInputFilepath();
        static void printHelpText();
        static bool file_exists(const std::__cxx11::basic_string<char> &name);
};


#endif //CONTEXTSV_CLI_H
