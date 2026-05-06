#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP

#include <string>

// Sostituisci "team_name" con il nome reale del tuo team
namespace team_name {

enum class AlignmentType {
    Global,
    Local,
    SemiGlobal
};

int Align(
    const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    AlignmentType type,
    int match,
    int mismatch,
    int gap,
    std::string* cigar = nullptr,
    unsigned int* target_begin = nullptr);

} // namespace team_name

#endif // ALIGNMENT_HPP