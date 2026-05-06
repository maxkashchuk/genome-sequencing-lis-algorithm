#ifndef MAIN
#define MAIN

#include <cstdint>
#include "bioparser/fasta_parser.hpp"

namespace Main {
    const std::string fasta_data_path_1 = "dataset/MAP006-1_2D_pass.fasta";
    
    const std::string fasta_data_path_2 = "dataset/MAP006-2_2D_pass.fasta";

    int opt;
    
    int match_score = -2;
    int mismatch_penalty = -1;
    int gap_penalty = -2;
    int k_mer_size = 3;
    int window_size = 3;

    std::string alignment_type = "global";

    struct Sequence {
        public:
            std::string name;
            std::string data;

            Sequence(
                const char* name_ptr, std::uint32_t name_len,
                const char* data_ptr, std::uint32_t data_len) : name(name_ptr, name_len),
                                                                data(data_ptr, data_len) {
                }
    };

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
                case 'k':
                    k_mer_size = std::atoi(optarg);
                    if (k_mer_size > 32) {
                        std::cerr << "Error. Maximum k-mer size is 32\n";
                        return -1;
                    }
                    break;
                case 'w':
                    window_size = std::atoi(optarg);
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
        std::cout << "K-mer size: " << k_mer_size << "\n";
        std::cout << "Window size: " << window_size << "\n";
    }

    void print_sequence(const Main::Sequence& seq) {
        std::cout << "Name: " << seq.name << "\n";
        std::cout << "Sequence: " << seq.data << "\n";
        std::cout << "Length of sequence: " << seq.data.size() << "\n";
    }

    // Srquence to binary template function (to be implemented)
    static void sequence_to_binary(std::vector<std::unique_ptr<Sequence>>& seq) {
        return;
    }
}

#endif