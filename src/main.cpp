//
// Created by perdomoj on 1/6/2023.
//

#include "cli.h"

#include <iostream>
#include <stdexcept>

int main(int argc, char *argv[]){
    // Run the command line interface
    CLI cli_obj;
    try
    {
        int exit_code;
        exit_code = cli_obj.parse(argc, argv);
        if (exit_code == 0)
        {
            // Run if arguments are provided
            cli_obj.run();
        }
    }

    catch (std::invalid_argument& e)
    {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    return 0;
}
