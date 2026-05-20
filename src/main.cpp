#include <iostream>
#include <string>
#include <vector>
#include <memory>
#include <unistd.h> 
#include <cstdlib>  

#include "main.h"
#include "minimizer.h"
#include "alignment.hpp"
#include "mapper.hpp"

int main(int argc, char* argv[]) {

    // 1. Parse command-line arguments
    uint32_t result = Main::arguments_parser(argc, argv);
    if (-1 == result) {
        return 0;
    }

    // 2. Load FASTA files into memory using bioparser
    std::vector<std::unique_ptr<Main::Sequence>> ref_seqs = 
        bioparser::Parser<Main::Sequence>::Create<bioparser::FastaParser>(Main::fasta_data_path_1)->Parse(-1);

    std::vector<std::unique_ptr<Main::Sequence>> query_seqs = 
        bioparser::Parser<Main::Sequence>::Create<bioparser::FastaParser>(Main::fasta_data_path_2)->Parse(-1);

    if (ref_seqs.empty() || query_seqs.empty()) {
        std::cerr << "Error: Unable to load sequences from the provided files.\n";
        return 1;
    }

    // 3. Map the alignment type
    team_name::AlignmentType type_enum;
    if (Main::alignment_type == "global") {
        type_enum = team_name::AlignmentType::Global;
    } else if (Main::alignment_type == "local") {
        type_enum = team_name::AlignmentType::Local;
    } else if (Main::alignment_type == "semi-global") {
        type_enum = team_name::AlignmentType::SemiGlobal;
    } else {
        std::cerr << "Error: Invalid alignment type. Use global, local, or semi-global.\n";
        return 1;
    }

    // 4. Initialize and extract minimizers for the reference genome
    MinimizerManager minimizer_mgr;
    std::cout << "[1/4] Extracting minimizers from the reference genome...\n";
    std::vector<Minimizer> ref_minimizers = minimizer_mgr.find_minimizers(ref_seqs, Main::k_mer_size, Main::window_size);

    // 5. Build the Mapper index
    Mapper mapper;
    std::cout << "[2/4] Building the hash index...\n";
    mapper.build_index(ref_minimizers, Main::frequency_threshold); // Filter out highly frequent k-mers

    // 6. Extract minimizers for the fragments (query)
    std::cout << "[3/4] Extracting minimizers from queries and searching for matches...\n";
    std::vector<Minimizer> query_minimizers = minimizer_mgr.find_minimizers(query_seqs, Main::k_mer_size, Main::window_size);

    // Identify matches and filter them using LIS
    std::vector<Match> raw_matches = mapper.find_matches(query_minimizers);
    std::vector<Match> lis_chain = mapper.filter_matches_lis(raw_matches);

    if (lis_chain.empty()) {
        std::cout << "No candidate region found between the sequences.\n";
        return 0;
    }

    // 7. Extract coordinates for precision alignment (Seed & Extend)
    std::cout << "[4/4] Executing precision alignment (Seed & Extend)...\n";

    // Retrieve the raw sequences (assuming for simplicity the first sequence in each file)
    std::string full_query = query_seqs[0]->data;
    std::string full_target = ref_seqs[0]->data;

    // The start and end coordinates are defined by the first and last match of the LIS chain.
    // Add k_mer_size to the end position to include the entire k-mer.
    int q_begin = lis_chain.front().query_pos;
    int q_end = lis_chain.back().query_pos + Main::k_mer_size;
    int t_begin = lis_chain.front().target_pos;
    int t_end = lis_chain.back().target_pos + Main::k_mer_size;

    // Safety boundary checks (Out of bounds)
    if (q_end > full_query.length()) q_end = full_query.length();
    if (t_end > full_target.length()) t_end = full_target.length();

    // Extract exclusively the candidate regions
    std::string query_region = full_query.substr(q_begin, q_end - q_begin);
    std::string target_region = full_target.substr(t_begin, t_end - t_begin);

    // 8. Execute dynamic programming alignment
    std::string cigar = "";
    unsigned int target_align_begin = 0;

    int score = team_name::Align(
        query_region.c_str(), query_region.length(),
        target_region.c_str(), target_region.length(),
        type_enum,
        Main::match_score,
        Main::mismatch_penalty,
        Main::gap_penalty,
        &cigar,
        &target_align_begin
    );

    // 9. Output results in a format inspired by PAF (Pairwise Mapping Format)
    std::cout << "\n--- Formatted PAF Output ---\n";
    std::cout << query_seqs[0]->name << "\t"        // 1. Query name
              << full_query.length() << "\t"        // 2. Query sequence length
              << q_begin << "\t"                    // 3. Query start
              << q_end << "\t"                      // 4. Query end
              << "+" << "\t"                        // 5. Relative strand
              << ref_seqs[0]->name << "\t"          // 6. Target name
              << full_target.length() << "\t"       // 7. Target sequence length
              << t_begin << "\t"                    // 8. Target start
              << t_end << "\t"                      // 9. Target end
              << score << "\t"                      // 10. Alignment Score (instead of matching bases)
              << std::max(query_region.length(), target_region.length()) << "\t" // 11. Alignment block length
              << "255\t"                            // 12. Mapping quality (255 = missing/uncalculated)
              << "cg:Z:" << cigar << "\n";          // 13. Optional CIGAR string

int main() {
    std::cout << "Hello, Bioinformatics!" << std::endl;
    return 0;
}