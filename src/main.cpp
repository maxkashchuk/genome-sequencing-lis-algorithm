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

    MinimizerManager minimizer_mgr;
    std::cout << "[1/4] Extracting minimizers from the reference genome...\n";
    
    std::vector<std::unique_ptr<Main::Sequence>> single_ref;
    single_ref.push_back(std::move(ref_seqs[0])); // Беремо тільки першу послідовність

    std::vector<Minimizer> ref_minimizers = minimizer_mgr.find_minimizers(single_ref, Main::k_mer_size, Main::window_size);

    // 5. Build the Mapper index
    Mapper mapper;
    std::cout << "[2/4] Building the hash index...\n";
    mapper.build_index(ref_minimizers, Main::frequency_threshold);

    std::cout << "[3/4] Processing queries one by one...\n";

    for (size_t i = 0; i < query_seqs.size(); ++i) {
        std::vector<std::unique_ptr<Main::Sequence>> current_query_vec;
        current_query_vec.push_back(std::move(query_seqs[i]));

        std::vector<Minimizer> query_minimizers = minimizer_mgr.find_minimizers(current_query_vec, Main::k_mer_size, Main::window_size);

        std::vector<Match> raw_matches = mapper.find_matches(query_minimizers);
        std::vector<Match> lis_chain = mapper.filter_matches_lis(raw_matches);

        if (lis_chain.empty()) {
            continue;
        }

        // --- Step 7 & 8 (Aligning reads) ---
        std::string full_query = current_query_vec[0]->data;
        std::string full_target = single_ref[0]->data;

        int q_begin = lis_chain.front().query_pos;
        int q_end = lis_chain.back().query_pos + Main::k_mer_size;
        int t_begin = lis_chain.front().target_pos;
        int t_end = lis_chain.back().target_pos + Main::k_mer_size;

        if (q_end > full_query.length()) q_end = full_query.length();
        if (t_end > full_target.length()) t_end = full_target.length();

        std::string query_region = full_query.substr(q_begin, q_end - q_begin);
        std::string target_region = full_target.substr(t_begin, t_end - t_begin);

        if (query_region.length() > 30000 || target_region.length() > 30000) {
            std::cerr << "Skipping read " << current_query_vec[0]->name << " - region too large for DP alignment\n";
            continue;
        }

        std::string cigar = "";
        unsigned int target_align_begin = 0;

        int score = team_name::Align(
            query_region.c_str(), query_region.length(),
            target_region.c_str(), target_region.length(),
            type_enum, Main::match_score, Main::mismatch_penalty, Main::gap_penalty,
            &cigar, &target_align_begin
        );

        // 9. Output PAF line for this specific read
        std::cout << current_query_vec[0]->name << "\t"        
                  << full_query.length() << "\t"        
                  << q_begin << "\t"                    
                  << q_end << "\t"                      
                  << "+" << "\t"                        
                  << single_ref[0]->name << "\t"
                  << full_target.length() << "\t"       
                  << t_begin << "\t"                    
                  << t_end << "\t"                      
                  << score << "\t"                      
                  << std::max(query_region.length(), target_region.length()) << "\t" 
                  << "255\t"                            
                  << "cg:Z:" << cigar << "\n";
    }

    return 0;
}