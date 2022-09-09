//
// Created by perdomoj on 9/8/2022.
//

#include <iostream>
#include <cstring>


int main(int argc, char *argv[]){
    // Print help text
    if(argc == 2 && (strcmp(argv[1], "-h")==0 || strcmp(argv[1], "--help")==0))
    {
        std::string help_text(
                "\n ===== GenSV =====\n"
                "\nPositional arguments:\n"
                "-i      BAM file input"
                "-o      Output directory\n"
                "\nOptional arguments:\n"
                "-h, --help     Show this help message and exit\n"
                "\nExample:\n"
                "./gensv -i path/to/input.bam -o /output_directory/");
        std::cout << help_text << std::endl;
    }
    return 0;
}
