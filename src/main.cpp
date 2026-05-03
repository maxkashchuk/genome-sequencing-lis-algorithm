#include <iostream>
#include <string>
#include <unistd.h> // Required for the getopt function
#include <cstdlib>  // Required for std::atoi (converts text to integers)

#include "main.h"
#include "minimizer.h"

#include <typeinfo>

int main(int argc, char* argv[]) {

    uint32_t result = Main::arguments_parser(argc, argv);
    if (-1 == result) {
        return 0;
    }

    std::vector<std::unique_ptr<Main::Sequence>> s = bioparser::Parser<Main::Sequence>::Create<bioparser::FastaParser>(Main::fasta_data_path)->Parse(-1);

    // Uncomment the line below to print the first sequence's name and data
    // Main::print_sequence(*s[0]);

    MinimizerManager minimizer;
    std::vector<Minimizer> minimizers = minimizer.find_minimizers(s, Main::k_mer_size, Main::window_size);

    std::cout << "Found " << minimizers.size() << " minimizers.\n";

    return 0;
}