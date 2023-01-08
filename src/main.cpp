//
// Created by perdomoj on 1/6/2023.
//

#include "CLI.h"

int main(int argc, char *argv[]){
    // Run the command line interface
    CLI cli_obj;
    cli_obj.parse(argc, argv);

    return 0;
}
