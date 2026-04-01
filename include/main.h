#ifndef MAIN
#define MAIN

#include <cstdint>
#include "bioparser/fasta_parser.hpp"

namespace Main {
    const std::string fasta_data_path = "dataset/MAP006-1_2D_pass.fasta";
    
    int opt;
    
    // Variables to store the values decoded by the parser
    int match_score = 0;
    int mismatch_penalty = 0;
    int gap_penalty = 0;
    std::string alignment_type = "";

    struct Sequence {
        public:
            Sequence(
                const char*, std::uint32_t,
                const char*, std::uint32_t) {
            }
    };
}

#endif