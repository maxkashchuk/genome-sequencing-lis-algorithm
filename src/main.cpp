#include <iostream>
#include <string>
#include <unistd.h> // Required for the getopt function
#include <cstdlib>  // Required for std::atoi (converts text to integers)

int main(int argc, char* argv[]) {
    int opt;
    
    // Variables to store the values decoded by the parser
    int match_score = 0;
    int mismatch_penalty = 0;
    int gap_penalty = 0;
    std::string alignment_type = "";

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
                return 1; // Terminates the program with an error code
        }
    }

    // At this point, the parsing is complete. 
    // You can pass these variables to the functions handling the bioinformatics logic.
    std::cout << "Parameters acquired successfully:\n";
    std::cout << "Type: " << alignment_type << " | Match: " << match_score 
              << " | Mismatch: " << mismatch_penalty << " | Gap: " << gap_penalty << "\n";

    return 0;
}