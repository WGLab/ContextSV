//
// Created by perdomoj on 1/6/2023.
//

#include "cli.h"

#include <iostream>
#include <stdexcept>

int main(int argc, char *argv[]){
    // Run the command line interface
    CLI cli;
    try
    {
        // Parse the input arguments and set the common parameters
        int exit_code;
        exit_code = cli.parse(argc, argv);
        if (exit_code == 0)
        {
            // Run the integrative caller
            cli.run();
        }
    }

    catch (std::invalid_argument& e)
    {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    return 0;
}
