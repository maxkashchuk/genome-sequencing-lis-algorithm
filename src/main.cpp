#include <iostream>
#include <string>
#include <unistd.h> // Required for the getopt function
#include <cstdlib>  // Required for std::atoi (converts text to integers)

#include "main.h"
#include "minimizer.h"
#include "alignment.hpp"

#include <typeinfo>

int main(int argc, char* argv[]) {

    uint32_t result = Main::arguments_parser(argc, argv);
    if (-1 == result) {
        return 0;
    }

    std::vector<std::unique_ptr<Main::Sequence>> s_1 = bioparser::Parser<Main::Sequence>::Create<bioparser::FastaParser>(Main::fasta_data_path_1)->Parse(-1);

    std::vector<std::unique_ptr<Main::Sequence>> s_2 = bioparser::Parser<Main::Sequence>::Create<bioparser::FastaParser>(Main::fasta_data_path_2)->Parse(-1);

    team_name::AlignmentType type_enum;
    if (Main::alignment_type == "global") {
        type_enum = team_name::AlignmentType::Global;
    } else if (Main::alignment_type == "local") {
        type_enum = team_name::AlignmentType::Local;
    } else if (Main::alignment_type == "semi-global") {
        type_enum = team_name::AlignmentType::SemiGlobal;
    } else {
        std::cerr << "Errore: Tipo non valido. Usa global, local, semi-global.\n";
        return 1;
    }

    std::string cigar = "";
    unsigned int target_begin = 0;

    std::string query = s_1[0]->data;
    std::string target = s_2[0]->data;

    int score = team_name::Align(
        query.c_str(), query.length(),
        target.c_str(), target.length(),
        type_enum,
        Main::match_score,
        Main::mismatch_penalty,
        Main::gap_penalty,
        &cigar,
        &target_begin
    );

    std::cout << "--- Alignment Results ---\n";
    std::cout << "Type       : " << Main::alignment_type << "\n";
    std::cout << "Score      : " << score << "\n";
    std::cout << "Target Pos : " << target_begin << "\n";
    std::cout << "CIGAR      : " << cigar << "\n";

    MinimizerManager minimizer_1;
    std::vector<Minimizer> minimizers_1 = minimizer_1.find_minimizers(s_1, Main::k_mer_size, Main::window_size);

    std::cout << "Found first sequence " << minimizers_1.size() << " minimizers.\n";

    MinimizerManager minimizer_2;
    std::vector<Minimizer> minimizers_2 = minimizer_2.find_minimizers(s_2, Main::k_mer_size, Main::window_size);

    std::cout << "Found second sequence " << minimizers_2.size() << " minimizers.\n";

    return 0;
}