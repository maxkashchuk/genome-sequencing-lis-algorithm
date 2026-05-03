#ifndef MINIMIZER
#define MINIMIZER

#include <vector>
#include <deque>
#include "main.h"

struct Minimizer {
    uint64_t value;
    int position;
};

class MinimizerManager {
    public:
        MinimizerManager() {}
    
    inline uint64_t encode(char c) {
        switch (c) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
        }
        return 0;
    }
    
    std::vector<Minimizer> find_minimizers(std::vector<std::unique_ptr<Main::Sequence>>& sequence, int k, int w) {
        std::vector<Minimizer> all_minimizers;

        uint64_t mask = (k == 32) ? ~0ULL : ((1ULL << (2 * k)) - 1);

        for (auto& seq_ptr : sequence) {
            auto& seq = *seq_ptr;
            const std::string& s = seq.data;

            int n = s.size();
            if (n < k) {
                continue;
            }

            std::vector<uint64_t> kmers;
            kmers.reserve(n - k + 1);

            uint64_t kmer = 0;

            for (int i = 0; i < n; i++) {
                kmer = ((kmer << 2) | encode(s[i])) & mask;

                if (i >= k - 1) {
                    kmers.push_back(kmer);
                }
            }

            std::deque<int> dq;

            for (int i = 0; i < kmers.size(); i++) {

                while (!dq.empty() && dq.front() <= i - w)
                    dq.pop_front();

                while (!dq.empty() && kmers[dq.back()] >= kmers[i])
                    dq.pop_back();

                dq.push_back(i);

                if (i >= w - 1) {
                    all_minimizers.push_back({kmers[dq.front()], dq.front()});
                }
            }
        }
        
        return all_minimizers;
    }
};

#endif