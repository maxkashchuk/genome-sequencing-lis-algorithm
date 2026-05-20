#ifndef MAPPER_H
#define MAPPER_H

#include <vector>
#include <unordered_map>
#include <algorithm>
#include "minimizer.h" 
#include "main.h"      

// Represents an exact match between a minimizer from the read (query) and the reference genome (target).
struct Match {
    int query_pos;
    int target_pos;

    // Defines the comparison operator to easily sort Match objects.
    // Primary sorting is by query_pos (ascending). If two matches share the same query_pos,
    // target_pos is used as a secondary criterion (descending).
    bool operator<(const Match& other) const {
        if (query_pos != other.query_pos) {
            return query_pos < other.query_pos; // Ascending order for the query
        }
        return target_pos > other.target_pos; // Descending order for the target
    }
};

class Mapper {
private:
    // Reference genome index:
    // Key: Minimizer hash value (uint64_t)
    // Value: List of positions (int) where the minimizer occurs in the reference genome
    std::unordered_map<uint64_t, std::vector<int>> reference_index;

public:
    Mapper() = default;

    // 1. Builds the index starting from the reference minimizers.
    // The frequency threshold parameter 'f' is used to identify and discard highly repetitive minimizers.
    void build_index(const std::vector<Minimizer>& ref_minimizers, double f = 0.001);

    // 2. Queries the index to find all raw matches for the given sequence fragment.
    std::vector<Match> find_matches(const std::vector<Minimizer>& query_minimizers);

    // 3. Core algorithm: applies the Longest Increasing Subsequence (LIS) in O(n log n) time.
    // Returns the filtered and sorted collinear chain representing the candidate region for alignment.
    std::vector<Match> filter_matches_lis(std::vector<Match>& matches);
};

#endif // MAPPER_H