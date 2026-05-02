#ifndef MAIN
#define MAIN

#include <cstdint>
#include "bioparser/fasta_parser.hpp"

namespace Main {
    const std::string fasta_data_path = "dataset/MAP006-1_2D_pass.fasta";
    
    int opt;
    
    int match_score = 0;
    int mismatch_penalty = 0;
    int gap_penalty = 0;
    std::string alignment_type = "";

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

    // Srquence to binary template function (to be implemented)
    void sequence_to_binary(std::vector<std::unique_ptr<Sequence>>& seq) {
        return;
    }
}

#endif