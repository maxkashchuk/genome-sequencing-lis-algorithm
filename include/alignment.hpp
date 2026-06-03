#ifndef ALIGNMENT_HPP
#define ALIGNMENT_HPP

#include <string>

namespace SequenceAlignment {

// Defines the algorithmic approach for the sequence alignment
enum class AlignmentType {
    Global,     // Aligns the entire sequence (e.g., Needleman-Wunsch)
    Local,      // Finds the best matching subsequences (e.g., Smith-Waterman)
    SemiGlobal  // Aligns the entire query to a portion of the target without end-gap penalties
};

// Performs the sequence alignment based on the specified parameters.
// Returns the optimal alignment score.
int Align(
    const char* query,                   // Pointer to the query sequence character array
    unsigned int query_len,              // Length of the query sequence
    const char* target,                  // Pointer to the target sequence character array
    unsigned int target_len,             // Length of the target sequence
    AlignmentType type,                  // The type of alignment to perform (Global, Local, or SemiGlobal)
    int match,                           // Positive score reward for matching characters
    int mismatch,                        // Penalty score for mismatching characters
    int gap,                             // Penalty score for introducing a gap
    std::string* cigar = nullptr,        // Optional pointer to store the resulting CIGAR string
    unsigned int* target_begin = nullptr // Optional pointer to store the 0-based start index on the target
);

} // namespace SequenceAlignment

#endif // ALIGNMENT_HPP
