#include "alignment.hpp"
#include <vector>
#include <algorithm>

namespace team_name {

int Align(
    const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    AlignmentType type,
    int match,
    int mismatch,
    int gap,
    std::string* cigar,
    unsigned int* target_begin)
{
    size_t n = query_len;
    size_t m = target_len;

    // Vettore 1D per simulare la matrice 2D (Cache-friendly)
    std::vector<int> H((n + 1) * (m + 1), 0);
    std::vector<char> trace((n + 1) * (m + 1), 0); // 'D'=Diag, 'U'=Up, 'L'=Left

    auto idx = [&](size_t i, size_t j) { return i * (m + 1) + j; };

    // Inizializzazione della matrice in base al tipo di allineamento
    for (size_t i = 1; i <= n; ++i) {
        if (type == AlignmentType::Global) {
            H[idx(i, 0)] = i * gap;
            trace[idx(i, 0)] = 'U';
        } else {
            H[idx(i, 0)] = 0;
        }
    }
    for (size_t j = 1; j <= m; ++j) {
        if (type == AlignmentType::Global) {
            H[idx(0, j)] = j * gap;
            trace[idx(0, j)] = 'L';
        } else {
            H[idx(0, j)] = 0;
        }
    }

    int max_score = -999999;
    size_t max_i = 0, max_j = 0;

    // Popolamento della matrice
    for (size_t i = 1; i <= n; ++i) {
        for (size_t j = 1; j <= m; ++j) {
            int score_diag = H[idx(i - 1, j - 1)] + (query[i - 1] == target[j - 1] ? match : mismatch);
            int score_up   = H[idx(i - 1, j)] + gap;
            int score_left = H[idx(i, j - 1)] + gap;

            int current_score = score_diag;
            char dir = 'D'; // Diagonale: M (Match/Mismatch)

            if (score_up > current_score) {
                current_score = score_up;
                dir = 'U'; // Su: I (Inserzione nel read)
            }
            if (score_left > current_score) {
                current_score = score_left;
                dir = 'L'; // Sinistra: D (Delezione dal read)
            }

            // Per l'allineamento Locale, il punteggio non può scendere sotto lo zero
            if (type == AlignmentType::Local && current_score < 0) {
                current_score = 0;
                dir = '0'; 
            }

            H[idx(i, j)] = current_score;
            trace[idx(i, j)] = dir;

            // Salvataggio del picco massimo per l'allineamento Locale
            if (type == AlignmentType::Local && current_score > max_score) {
                max_score = current_score;
                max_i = i;
                max_j = j;
            }
        }
    }

    // Ricerca del punto di inizio per il Backtracking
    size_t start_i = n, start_j = m;
    int final_score = 0;

    if (type == AlignmentType::Global) {
        final_score = H[idx(n, m)];
    } else if (type == AlignmentType::Local) {
        final_score = max_score;
        start_i = max_i;
        start_j = max_j;
    } else if (type == AlignmentType::SemiGlobal) {
        final_score = H[idx(n, 0)];
        start_i = n; start_j = 0;
        for (size_t j = 1; j <= m; ++j) {
            if (H[idx(n, j)] > final_score) {
                final_score = H[idx(n, j)];
                start_i = n; start_j = j;
            }
        }
        for (size_t i = 1; i <= n; ++i) {
            if (H[idx(i, m)] > final_score) {
                final_score = H[idx(i, m)];
                start_i = i; start_j = m;
            }
        }
    }

    // Backtracking per calcolare la CIGAR e la posizione di inizio
    size_t curr_i = start_i;
    size_t curr_j = start_j;
    std::string raw_cigar = "";

    while (curr_i > 0 || curr_j > 0) {
        if (type == AlignmentType::Local && H[idx(curr_i, curr_j)] == 0) break;
        if (type == AlignmentType::SemiGlobal && (curr_i == 0 || curr_j == 0)) break;

        char dir = trace[idx(curr_i, curr_j)];
        if (dir == 'D') {
            raw_cigar += 'M';
            curr_i--;
            curr_j--;
        } else if (dir == 'U') {
            raw_cigar += 'I'; 
            curr_i--;
        } else if (dir == 'L') {
            raw_cigar += 'D'; 
            curr_j--;
        } else {
            break;
        }
    }

    // Assegnazione dei parametri opzionali richiesti dalla consegna
    if (target_begin != nullptr) {
        *target_begin = curr_j; 
    }

    if (cigar != nullptr) {
        std::reverse(raw_cigar.begin(), raw_cigar.end());
        std::string compact_cigar = "";
        if (!raw_cigar.empty()) {
            int count = 1;
            for (size_t i = 1; i < raw_cigar.length(); ++i) {
                if (raw_cigar[i] == raw_cigar[i - 1]) {
                    count++;
                } else {
                    compact_cigar += std::to_string(count) + raw_cigar[i - 1];
                    count = 1;
                }
            }
            compact_cigar += std::to_string(count) + raw_cigar.back();
        }
        *cigar = compact_cigar;
    }

    return final_score;
}

} // namespace team_name