#include <iostream>
#include <string>
#include <unistd.h> // Required for the getopt function
#include <cstdlib>  // Required for std::atoi (converts text to integers)

#include "main.h"

uint32_t arguments_parser(int argc, char* argv[])
{
    using namespace Main;

     // The third parameter "a:m:n:g:" defines which flags to look for.
    // The colon (:) indicates that the preceding flag requires an associated value (e.g., -m 5)
    
    while ((opt = getopt(argc, argv, "a:m:n:g:")) != -1) {
        switch (opt) {
            case 'a':
                alignment_type = optarg; // optarg is a global variable containing the text following the flag
                break;
            case 'm':
                match_score = std::atoi(optarg); // Converts the text string into an integer
                break;
            case 'n':
                mismatch_penalty = std::atoi(optarg);
                break;
            case 'g':
                gap_penalty = std::atoi(optarg);
                break;
            default:
                std::cerr << "Error. Correct usage: " << argv[0] << " -a <type> -m <match> -n <mismatch> -g <gap>\n";
                return -1; // Terminates the program with an error code
        }
    }

    return 1; // Indicates successful parsing of arguments
}

void print_arguments() {
    using namespace Main;
    std::cout << "Alignment Type: " << alignment_type << "\n";
    std::cout << "Match Score: " << match_score << "\n";
    std::cout << "Mismatch Penalty: " << mismatch_penalty << "\n";
    std::cout << "Gap Penalty: " << gap_penalty << "\n";
}

int main(int argc, char* argv[]) {

    uint32_t result = arguments_parser(argc, argv);
    if (-1 == result) {
        return 0;
    }

    return 0;
}