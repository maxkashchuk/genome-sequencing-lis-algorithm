#ifndef MAIN
#define MAIN

#include <cstdint>
#include <string>
#include <vector>
#include <memory>
#include <iostream>
#include <unistd.h> 
#include <cstdlib>  
#include "bioparser/fasta_parser.hpp"

namespace Main {
    const std::string reference_data_path_1 = "dataset/GCF_000005845.2_ASM584v2_genomic.fna";
    const std::string query_data_path_2 = "dataset/MAP006-1_2D_pass.fasta";

    inline int opt;
    
    // Alignment parameters
    inline int match_score = -2;
    inline int mismatch_penalty = -1;
    inline int gap_penalty = -2;
    inline std::string alignment_type = "global";

    // Minimizer parameters with standard default values
    inline int k_mer_size = 15;
    inline int window_size = 5;
    inline double frequency_threshold = 0.001;

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

    inline int arguments_parser(int argc, char* argv[])
    {
        // The string "a:m:n:g:k:w:f:" specifies to getopt which flags require an associated argument
        while ((opt = getopt(argc, argv, "a:m:n:g:k:w:f:")) != -1) {
            switch (opt) {
                case 'a':
                    alignment_type = optarg;
                    break;
                case 'm':
                    match_score = std::atoi(optarg);
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
                        std::cerr << "Error. The maximum allowed k-mer size is 32.\n";
                        return -1;
                    }
                    break;
                case 'w':
                    window_size = std::atoi(optarg);
                    break;
                case 'f':
                    frequency_threshold = std::atof(optarg);
                    if (frequency_threshold < 0.0 || frequency_threshold > 1.0) {
                        std::cerr << "Error. The frequency threshold must be between 0.0 and 1.0.\n";
                        return -1;
                    }
                    break;
                default:
                    std::cerr << "Usage: " << argv[0] 
                              << " -a <type> -m <match> -n <mismatch> -g <gap> -k <kmer> -w <window> -f <freq>\n";
                    return -1;
            }
        }
        return 1;
    }

    inline void print_arguments() {
        std::cout << "--- Configuration Parameters ---\n";
        std::cout << "Alignment Type   : " << alignment_type << "\n";
        std::cout << "Match Score      : " << match_score << "\n";
        std::cout << "Mismatch Penalty : " << mismatch_penalty << "\n";
        std::cout << "Gap Penalty      : " << gap_penalty << "\n";
        std::cout << "K-mer Size       : " << k_mer_size << "\n";
        std::cout << "Window Size      : " << window_size << "\n";
        std::cout << "Frequency Thresh.: " << frequency_threshold << "\n";
        std::cout << "--------------------------------\n";
    }

    inline void print_sequence(const Main::Sequence& seq) {
        std::cout << "Name: " << seq.name << "\n";
        std::cout << "Sequence length: " << seq.data.size() << "\n";
    }
}

#endif // MAIN