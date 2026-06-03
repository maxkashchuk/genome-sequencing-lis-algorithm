#include "mapper.hpp"
#include <iostream>
#include <algorithm>

// Build the reference index
void Mapper::build_index(const std::vector<Minimizer>& ref_minimizers, double f) {
    // Calculate frequencies to filter out overly abundant minimizers
    std::unordered_map<uint64_t, int> counts;
    for (const auto& m : ref_minimizers) {
        counts[m.value]++;
    }

    int total_minimizers = ref_minimizers.size();
    for (const auto& m : ref_minimizers) {
        // Ignore the minimizer if its occurrence frequency exceeds the threshold 'f'
        if (counts[m.value] > (total_minimizers * f)) continue;
        reference_index[m.value].push_back(m.position);
    }
}

// Find matches between the query and the reference
std::vector<Match> Mapper::find_matches(const std::vector<Minimizer>& query_minimizers) {
    std::vector<Match> matches;
    for (const auto& qm : query_minimizers) {
        if (reference_index.count(qm.value)) {
            for (int ref_pos : reference_index[qm.value]) {
                matches.push_back({qm.position, ref_pos});
            }
        }
    }
    // Sort matches by their position on the query to prepare for the LIS algorithm
    std::sort(matches.begin(), matches.end());
    return matches;
}

// LIS (Longest Increasing Subsequence) Algorithm in O(n log n)
std::vector<Match> Mapper::filter_matches_lis(std::vector<Match>& matches) {
    if (matches.empty()) return {};

    size_t n = matches.size();
    std::vector<int> tails; // Array storing the optimal target positions (LI)
    std::vector<int> parent(n, -1);
    std::vector<int> tails_indices; // Indices of the matches stored in 'tails'

    for (int i = 0; i < n; ++i) {
        // Find the correct position to insert or update the current match using binary search
        auto it = std::lower_bound(tails.begin(), tails.end(), matches[i].target_pos);
        int pos = std::distance(tails.begin(), it);

        if (it == tails.end()) {
            tails.push_back(matches[i].target_pos);
            tails_indices.push_back(i);
        } else {
            *it = matches[i].target_pos;
            tails_indices[pos] = i;
        }

        // Store the parent index to enable backtracking
        if (pos > 0) {
            parent[i] = tails_indices[pos - 1];
        }
    }

    // Reconstruct the optimal sequence chain via backtracking
    std::vector<Match> result;
    int curr = tails_indices.back();
    while (curr != -1) {
        result.push_back(matches[curr]);
        curr = parent[curr];
    }
    std::reverse(result.begin(), result.end());
    
    return result;
}
