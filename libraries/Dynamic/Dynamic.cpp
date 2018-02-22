//
// Created by matep on 03/12/2017.
//
#include "Dynamic.h"
#include "Common.hpp"
#include <iostream>
#include <iomanip>
#include <cstring>
#include <algorithm>

using namespace std;


// Binary search (note boundaries in the caller)
int cell_index(std::vector<int> &v, int l, int r, int key) {
    while (r - l > 1) {
        int m = l + (r - l) / 2;
        if (v[m] >= key)
            r = m;
        else
            l = m;
    }

    return r;
}

/**
 * Finds the longest increasing subsequence in a vector of integers.
 * @param v A reference to a vector of integers.
 * @return The number of elements in the longest increasing subsequence.
 */
int lis(const std::vector<int> &v) {
    int v_size = v.size();
    if (v_size == 0)
        return 0;

    std::vector<int> tail(v_size, 0);
    int length = 1; // always points empty slot in tail

    tail[0] = v[0];
    for (size_t i = 1; i < v_size; i++) {
        if (v[i] < tail[0])
            // new smallest value
            tail[0] = v[i];
        else if (v[i] > tail[length - 1])
            // v[i] extends largest subsequence
            tail[length++] = v[i];
        else
            // v[i] will become end candidate of an existing subsequence or
            // Throw away larger elements in all LIS, to make room for upcoming grater elements than v[i]
            // (and also, v[i] would have already appeared in one of LIS, identify the location and replace it)
            tail[cell_index(tail, -1, length - 1, v[i])] = v[i];
    }

    return length;
}


/**
 * Finds all overlaps of a given query sequence based on the longest increasing subsequence algorithm and a specified threshold.
 * @param query_id The id of the query sequence.
 * @param minimizers The query sequence's minimizers.
 * @param lookup_map A map providing all minimizers (collected from all seqeunces) for a given hash.
 * @param lis_threshold The LIS threshold which must be satisfied in order to classify two sequences as overlaps
 * @param occurrences A map mapping a minimizer hash to its number of occurences (in all sequences)
 * @return All overlaps of a given query seqeunce, denoted by the IDs of the found overlapping sequences as well as a boolean specifying the strand (true for same strand, false for different)
 */
vector<pair<int, bool>> find_overlaps_by_LIS(int query_id,
                                             vector<minim> &minimizers,
                                             unordered_map<uint64_t, vector<hashMinPair3>> &lookup_map,
                                             int lis_threshold,
                                             unordered_map<uint64_t, uint32_t> &occurrences) {
    unordered_map<uint64_t, vector<int>> same_strand;
    unordered_map<uint64_t, vector<int>> different_strand;

    for (auto min : minimizers) {
        if (occurrences[min.hash] > IGNORE_THRESHOLD)
            continue;

        auto matches = lookup_map.find(min.hash);
        if (matches == lookup_map.end())
            continue;

        for (auto match : matches->second) {
            if (match.seq_id > query_id) {
                continue;
            }

            int final = abs(match.index);
            if (match.index * min.index < 0) {
                different_strand[match.seq_id].push_back(-final);
            } else {
                same_strand[match.seq_id].push_back(final);
            }
        }
    }

    vector<pair<int, bool>> overlaps;
    for (auto &entry : same_strand) {
        if (lis(entry.second) >= lis_threshold) {
            overlaps.emplace_back(entry.first, true);
        }
    }
    for (auto &entry : different_strand) {
        if (lis(entry.second) >= lis_threshold) {
            overlaps.emplace_back(make_pair(entry.first, false));
        }
    }

    return overlaps;
};
